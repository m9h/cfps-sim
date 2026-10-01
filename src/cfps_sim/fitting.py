"""Bounded log-space least squares with GPU residuals and linear algebra."""
from types import SimpleNamespace
import torch


def least_squares(residual, start, bounds, max_iterations=100):
    """Projected damped Gauss--Newton with bounded central differences.

    Numerical derivatives are checked at two step sizes. Model
    evaluations, residuals, Jacobians, and linear solves stay on the GPU.
    """
    lower, upper = bounds
    x = start.clone().clamp(lower, upper)
    damping = 1e-3
    success = False
    termination = "iteration_limit"
    evaluations = 0

    def evaluate(value):
        nonlocal evaluations
        evaluations += 1
        result = residual(value)
        if not torch.isfinite(result).all():
            raise RuntimeError("Nonfinite calibration residuals")
        return result

    def jacobian(value, step=1e-4):
        columns = []
        for index in range(len(value)):
            plus, minus = value.clone(), value.clone()
            plus[index] = (plus[index] + step).clamp(lower[index], upper[index])
            minus[index] = (minus[index] - step).clamp(lower[index], upper[index])
            columns.append((evaluate(plus)-evaluate(minus))/(plus[index]-minus[index]))
        return torch.stack(columns, dim=1)

    r = evaluate(x)
    cost = r.square().sum()/2
    for iteration in range(max_iterations):
        jac = jacobian(x)
        gradient = jac.T @ r
        projected = gradient.clone()
        projected[(x <= lower) & (gradient > 0)] = 0
        projected[(x >= upper) & (gradient < 0)] = 0
        normal = jac.T @ jac
        # Scale optimality by channel sensitivity and residual magnitude. An
        # absolute gradient threshold depends on the supplied measurement SDs.
        scaled_gradient = (projected / normal.diagonal().clamp_min(1e-30).sqrt()
                           / r.norm().clamp_min(1)).abs().max()
        if scaled_gradient < 1e-6:
            success = True
            termination = "scaled_gradient"
            break
        scale = normal.diagonal().clamp_min(1e-12)
        accepted = False
        for attempt in range(20):
            delta = torch.linalg.solve(normal + damping*torch.diag(scale), -gradient)
            candidate = (x + delta).clamp(lower, upper)
            trial = evaluate(candidate)
            trial_cost = trial.square().sum()/2
            if trial_cost < cost:
                improvement = cost-trial_cost
                step_size = (candidate-x).abs().max()
                actual_step = candidate-x
                predicted_improvement = -gradient@actual_step - actual_step@(normal@actual_step)/2
                x, r, cost = candidate, trial, trial_cost
                damping = max(damping/3, 1e-12)
                accepted = True
                if improvement < 1e-9*(1+cost) and abs(predicted_improvement) < 1e-9*(1+cost):
                    success = True
                    termination = "objective_tolerance"
                elif step_size < 1e-8*(1+x.abs().max()):
                    success = True
                    termination = "step_tolerance"
                break
            damping *= 10
        if success or not accepted:
            if not success:
                termination = "no_descent_step"
            break
    jac = jacobian(x)
    refined = jacobian(x, step=5e-5)
    derivative_difference = (jac-refined).norm()/refined.norm().clamp_min(1e-30)
    if derivative_difference > 1e-3:
        success = False
        termination = "unstable_jacobian"
    return SimpleNamespace(x=x, jac=refined, cost=cost, success=success,
                           termination=termination, scaled_gradient=float(scaled_gradient),
                           iterations=iteration+1, evaluations=evaluations,
                           derivative_relative_difference=float(derivative_difference))
