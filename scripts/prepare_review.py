"""Build a review packet that exposes blocked gates without proposing wet-lab runs."""
import json
from pathlib import Path
from cfps_sim.preflight import audit_recipe

ROOT = Path(__file__).resolve().parents[1]

if __name__ == '__main__':
    registry = json.loads((ROOT/'data/registry/inventory.draft.json').read_text())
    recipe = json.loads((ROOT/'examples/recipe.json').read_text())
    prices = json.loads((ROOT/'examples/prices.draft.json').read_text())
    audit = audit_recipe(recipe, registry, prices)
    packet = {
        'status': 'pre_experiment_review_not_submission',
        'registry_stocks': len(registry['stocks']),
        'price_conflicts': registry['price_conflicts'],
        'selected_accounting_source': registry.get('selected_accounting_source'),
        'draft_ambiguities': registry.get('draft_ambiguities', []),
        'local_calibration_plan': json.loads((ROOT/'experiments/local-calibration.json').read_text()),
        'illustrative_recipe_audit': audit,
        'prioritized_hypotheses': [
            {'priority': 1, 'module': 'salt_buffer_nucleotide_balance',
             'question': 'Does Mg/NTP balance or buffering explain failed formulations?',
             'missing_evidence': ['base 18-AA mixture reconciliation', 'available pH measurements', 'Mg/NTP response data'],
             'model_gap': 'analytical totals available; free Mg and pH effects absent'},
            {'priority': 2, 'module': 'energy_regeneration',
             'question': 'Is a formulation limited by energy substrate chemistry or effective enzyme capacity?',
             'missing_evidence': ['ATP/ADP or other energy readouts if available', 'stock activities', 'geometry-dependent transfer measurements'],
             'model_gap': 'lumped Fuel does not identify actual substrate chemistry'},
            {'priority': 3, 'module': 'folding_redox_activity',
             'question': 'Do PETase/reteplase readouts change through expression or active fraction?',
             'missing_evidence': ['target-specific activity calibration', 'orthogonal abundance measurement if available'],
             'model_gap': 'no validated target activity model'}],
        'shortlist_acceptance_gates': [
            'inventory, concentrations, execution constraints and price ledger reconciled',
            'all proposed stocks and mixed-stock constituents accounted for',
            'target construct version and assay conversion documented',
            'in-domain reference observations available, or panel explicitly designated calibration',
            'prediction performance compared with empirical baseline on complete held-out recipes',
            'candidate ranking stress-tested against fitted uncertainty and structural model alternatives',
            'plate layout includes references, controls, replication and batch metadata',
            'measurement requests and interpretation criteria specified before execution'],
        'candidate_rankings': None,
        'reason_no_rankings': 'The supplied draft completes this example price ledger, but in-domain calibration, operative rules and activity scoring remain unresolved; literature prediction is not a Ginkgo recipe score.'}
    (ROOT/'reports/pre_experiment_review.json').write_text(json.dumps(packet, indent=2)+'\n')
    print(f"Audited {packet['registry_stocks']} stocks; {len(packet['price_conflicts'])} price conflicts; shortlist blocked")
