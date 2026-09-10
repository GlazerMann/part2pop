from part2pop.aerosol_particle import Particle


def test_get_variable_forwards_keyword_arguments_for_critical_supersaturation(
    monkeypatch,
):
    particle = Particle.__new__(Particle)
    calls = []

    def fake_get_critical_supersaturation(*, T, return_D_crit=False):
        calls.append((T, return_D_crit))
        return (0.5, 1.2e-6) if return_D_crit else 0.5

    monkeypatch.setattr(
        particle,
        "get_critical_supersaturation",
        fake_get_critical_supersaturation,
    )

    result = particle.get_variable(
        "s_c",
        T=298.15,
        return_D_crit=True,
    )

    assert result == (0.5, 1.2e-6)
    assert calls == [(298.15, True)]
