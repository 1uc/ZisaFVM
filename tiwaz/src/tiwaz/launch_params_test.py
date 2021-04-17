import pytest

from tiwaz.launch_params import LaunchParams


def test_booleaness():
    lp = LaunchParams()
    assert not lp

    lp = LaunchParams([])
    assert not lp

    lp = LaunchParams([{"bla": "flu"}])
    assert lp


def test_iterable():
    values = ["a", "b", "c", "d"]
    lp = LaunchParams(values)
    assert lp == values


def test_product():
    a = [{"a": "1.1"}, {"a": "1.2"}]
    b = [{"b": "2.1"}, {"b": "2.2"}, {"b": "2.3"}]
    ab = [
        {"a": "1.1", "b": "2.1"},
        {"a": "1.1", "b": "2.2"},
        {"a": "1.1", "b": "2.3"},
        {"a": "1.2", "b": "2.1"},
        {"a": "1.2", "b": "2.2"},
        {"a": "1.2", "b": "2.3"},
    ]

    lp_a = LaunchParams(a)
    lp_b = LaunchParams(b)
    lp_ab = lp_a.product(lp_b)

    assert lp_ab == ab


def test_equality():
    ab = [
        {"a": "1.1", "b": "2.1"},
        {"a": "1.2", "b": "2.1"},
        {"a": "1.2", "b": "2.2"},
        {"a": "1.2", "b": "2.3"},
    ]
    ab_prime = [
        {"a": "1.1", "b": "2.1"},
        {"a": "1.2", "b": "2.2"},
        {"a": "1.2", "b": "2.1"},
        {"a": "1.2", "b": "2.3"},
    ]

    assert LaunchParams(ab) == LaunchParams(ab_prime)
    assert LaunchParams(ab) == ab_prime


def test_inequality():
    ab = [
        {"a": "1.1", "b": "2.1"},
        {"a": "1.2", "b": "2.1"},
        {"a": "1.2", "b": "2.4"},
        {"a": "1.2", "b": "2.3"},
    ]
    ab_prime = [
        {"a": "1.1", "b": "2.1"},
        {"a": "1.2", "b": "2.2"},
        {"a": "1.3", "b": "2.1"},
        {"a": "1.2", "b": "2.3"},
    ]

    assert LaunchParams(ab) != LaunchParams(ab_prime)
    assert LaunchParams(ab) != ab_prime


def test_lazy_product():
    a = [{"a": "1.1"}, {"a": "1.2"}]

    lp_none = LaunchParams()
    lp_a = LaunchParams(a)

    assert lp_none.product(lp_a) == lp_a
    assert lp_a.product(lp_none) == lp_a
