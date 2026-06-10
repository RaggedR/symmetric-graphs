"""Tests for pure-Python layout functions (no GAP dependency)."""
import sys, os, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))


def test_no_eval_in_gap_output():
    """Ensure draw_graph.py and draw_quotients.py use ast.literal_eval, not eval."""
    import inspect
    import draw_graph as dg
    import draw_quotients as dq
    src_dg = inspect.getsource(dg.run_gap)
    src_dq = inspect.getsource(dq.run_gap)
    assert 'eval(' not in src_dg or 'literal_eval(' in src_dg, \
        "eval() found in draw_graph.run_gap — use ast.literal_eval"
    assert 'eval(' not in src_dq or 'literal_eval(' in src_dq, \
        "eval() found in draw_quotients.run_gap — use ast.literal_eval"


def test_segments_cross_basic():
    """X-crossing segments should return True."""
    from draw_graph import segments_cross
    assert segments_cross((0,0),(1,1),(0,1),(1,0)) is True


def test_segments_no_cross_parallel():
    """Parallel non-overlapping segments should return False."""
    from draw_graph import segments_cross
    assert segments_cross((0,0),(1,0),(0,1),(1,1)) is False


def test_segments_no_cross_collinear():
    """Collinear non-overlapping segments should return False."""
    from draw_graph import segments_cross
    assert segments_cross((0,0),(1,0),(2,0),(3,0)) is False


def test_count_all_crossings_no_cross():
    """A simple square graph (cycle C4) with no crossings."""
    from draw_graph import count_all_crossings
    # Square layout: 4 vertices at unit square corners
    coords = {1: (0,0), 2: (1,0), 3: (1,1), 4: (0,1)}
    # adj[v-1] gives neighbors of vertex v (GAP convention: 0-indexed list, 1-indexed vertices)
    adj = [[2,4], [1,3], [2,4], [3,1]]
    assert count_all_crossings(coords, adj, 4, perturb=False) == 0
