from src.generate_panel import common_member

def test_common_member():
    assert common_member(["a", "b"], ["b", "c"]) == 1