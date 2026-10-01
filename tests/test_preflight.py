import json
from pathlib import Path

import pytest

from cfps_sim.preflight import audit_recipe, model_initial_conditions
from cfps_sim.benchmark import partition

ROOT = Path(__file__).resolve().parents[1]


def test_counterions_base_and_conflicts():
    recipe = json.loads((ROOT/'examples/recipe.json').read_text())
    registry = json.loads((ROOT/'data/registry/inventory.draft.json').read_text())
    prices = json.loads((ROOT/'examples/prices.draft.json').read_text())
    report = audit_recipe(recipe, registry, prices, {'magnesium': 2.6, 'potassium': 150})
    assert report['known_final_mM']['magnesium'] == pytest.approx(8.85)
    assert report['known_final_mM']['potassium'] == pytest.approx(237.5)
    assert report['known_final_mM']['glutamate'] == pytest.approx(100)
    assert report['known_final_mM']['atp'] == pytest.approx(1.5)
    assert not report['submission_ready']
    assert 'conflicting_price_sources' in report['blockers']
    assert report['cost']['total_cost_usd'] is None
    assert 'magnesium_glutamate' in report['unsupported_model_effects']
    with pytest.raises(ValueError, match='backgrounds'):
        model_initial_conditions(report, {}, {})


def test_unknown_base_and_unregistered_stocks_block():
    recipe = json.loads((ROOT/'examples/recipe.json').read_text())
    recipe['volumes_nl']['unknown_stock'] = 25
    recipe['volumes_nl']['nuclease_free_water'] -= 25
    registry = json.loads((ROOT/'data/registry/inventory.draft.json').read_text())
    report = audit_recipe(recipe, registry)
    assert 'base_buffer_composition_unknown' in report['blockers']
    assert report['absent_stocks'] == ['unknown_stock']
    with pytest.raises(ValueError, match='base composition'):
        model_initial_conditions(report, {}, {})


def test_recipe_replicates_cannot_leak_between_plate_blocks():
    rows = [{'plate': p, 'recipe_key': (p,)} for p in [1, 8, 9, 10]]
    rows.extend([{'plate': 1, 'recipe_key': (10,)}, {'plate': 8, 'recipe_key': (9,)}])
    blocks, purged = partition(rows)
    assert purged['train'] == 1
    assert purged['validation'] == 1
    sets = [set(row['recipe_key'] for row in block) for block in blocks.values()]
    for i, keys in enumerate(sets):
        assert all(not keys.intersection(other) for other in sets[i+1:])


def test_mapping_requires_complete_explicit_background_and_other_states():
    from cfps_sim.construct import AMINO_ACIDS
    recipe = json.loads((ROOT/'examples/recipe.json').read_text())
    registry = json.loads((ROOT/'data/registry/inventory.draft.json').read_text())
    audit = audit_recipe(recipe, registry, base_final_mM={})
    names = ['alanine', 'cysteine', 'aspartate', 'glutamate', 'phenylalanine', 'glycine',
             'histidine', 'isoleucine', 'lysine', 'leucine', 'methionine', 'asparagine',
             'proline', 'glutamine', 'arginine', 'serine', 'threonine', 'valine', 'tryptophan', 'tyrosine']
    backgrounds = {name: 0.1 for name in ['atp', 'gtp', 'ctp', 'utp'] + names}
    other = dict(ADP=0, GDP=0, Fuel=0, mRNA=0, Protein=0, GFP=0, Machinery=1)
    mapped = model_initial_conditions(audit, backgrounds, other)
    assert mapped['initial']['ATP'] == pytest.approx(1.6)
    assert mapped['initial']['aa_E'] == pytest.approx(100.1)
    assert {f'aa_{aa}' for aa in AMINO_ACIDS}.issubset(mapped['initial'])
    assert mapped['initial']['Fuel'] == 0
    assert mapped['recipe_blockers']
    with pytest.raises(ValueError, match='other_initial'):
        model_initial_conditions(audit, backgrounds, {})
