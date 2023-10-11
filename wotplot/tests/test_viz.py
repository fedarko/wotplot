# Tests the visualization functions. Note that testing matplotlib images is
# kind of a pain (https://stackoverflow.com/q/27948126), so these tests mostly
# boil down to checking that the correct logging output gets written -- this is
# sufficient for checking that we are going to the desired branches of the
# code, etc. There's more work that could be done to polish these tests, of
# course, but I think this should be sufficient for now.
import pytest
from matplotlib import pyplot
from wotplot import DotPlotMatrix, FWD, REV, BOTH
from wotplot._viz import viz_spy, viz_imshow


def test_viz_spy_notbinary_bad_draw_order():
    # same example as test_make.test_make_palindrome_not_binary()
    dpm = DotPlotMatrix("AATCGATC", "TATCGAT", 6, binary=False)
    exp_err = (
        f"draw_order must include exactly 3 elements ({FWD}, {REV}, and "
        f"{BOTH} in any order)."
    )
    bad_inputs = [
        (BOTH, BOTH, BOTH),
        (FWD,),
        (FWD, REV, BOTH, 3),
        [],
        (REV, FWD),
    ]
    for bi in bad_inputs:
        with pytest.raises(ValueError) as ei:
            viz_spy(dpm, draw_order=bi)
        assert str(ei.value) == exp_err


def _check_logging_output(exp_lines, out_lines):
    assert len(out_lines) == len(exp_lines)
    for ol, el in zip(out_lines, exp_lines):
        # this will break if any of the logged stuff after the time includes
        # "s:" but we're safe for now
        after_time = ol.split("s: ")[1]
        assert after_time == el


def test_viz_spy_notbinary_verbose(capsys):
    # It's hard to actually check the matplotlib graphical output, so for the
    # time being we just check here that the logging is done correctly.
    dpm = DotPlotMatrix("AATCGATC", "TATCGAT", 6, binary=False)
    viz_spy(dpm, verbose=True)
    exp_lines = (
        f'Visualizing "{FWD}" cells with spy()...\n'
        f'Done visualizing "{FWD}" cells.\n'
        f'Visualizing "{REV}" cells with spy()...\n'
        f'Done visualizing "{REV}" cells.\n'
        f'Visualizing "{BOTH}" cells with spy()...\n'
        f'Done visualizing "{BOTH}" cells.\n'
        f"Slightly restyling the visualization...\n"
        f"Done.\n"
    ).splitlines()
    out_lines = capsys.readouterr().out.splitlines()
    _check_logging_output(exp_lines, out_lines)


def test_viz_spy_notbinary_verbose_diff_bgcolor(capsys):
    dpm = DotPlotMatrix("AATCGATC", "TATCGAT", 6, binary=False)
    # this is actually the worst color scheme known to man, but it works for
    # this test (i wanna make sure that my code recognizes that the requested
    # bg color isn't white and accordingly updates it)
    viz_spy(
        dpm, nbcmap={0: "#000", 1: "#ff0", -1: "#f0f", 2: "#fff"}, verbose=True
    )
    exp_lines = (
        f"Setting background color to #000...\n"
        f"Done setting background color.\n"
        f'Visualizing "{FWD}" cells with spy()...\n'
        f'Done visualizing "{FWD}" cells.\n'
        f'Visualizing "{REV}" cells with spy()...\n'
        f'Done visualizing "{REV}" cells.\n'
        f'Visualizing "{BOTH}" cells with spy()...\n'
        f'Done visualizing "{BOTH}" cells.\n'
        f"Slightly restyling the visualization...\n"
        f"Done.\n"
    ).splitlines()
    out_lines = capsys.readouterr().out.splitlines()
    _check_logging_output(exp_lines, out_lines)


def test_viz_spy_binary_verbose(capsys):
    dpm = DotPlotMatrix("AATCGATC", "TATCGAT", 6)
    viz_spy(dpm, verbose=True)
    exp_lines = (
        "Visualizing all match cells with spy()...\n"
        "Done visualizing all match cells.\n"
        "Slightly restyling the visualization...\n"
        "Done.\n"
    ).splitlines()
    out_lines = capsys.readouterr().out.splitlines()
    _check_logging_output(exp_lines, out_lines)


def test_viz_imshow_binary_verbose(capsys):
    dpm = DotPlotMatrix("AATCGATC", "TATCGAT", 6)
    viz_imshow(dpm, verbose=True)
    exp_lines = (
        "Converting the matrix to dense format...\n"
        "Calling imshow()...\n"
        "Slightly restyling the visualization...\n"
        "Done.\n"
    ).splitlines()
    out_lines = capsys.readouterr().out.splitlines()
    _check_logging_output(exp_lines, out_lines)


def test_viz_imshow_notbinary_verbose(capsys):
    dpm = DotPlotMatrix("AATCGATC", "TATCGAT", 6, binary=False)
    viz_imshow(dpm, verbose=True)
    exp_lines = (
        "Converting the matrix to dense format...\n"
        "Converting the matrix from numbers to colors...\n"
        "Calling imshow()...\n"
        "Slightly restyling the visualization...\n"
        "Done.\n"
    ).splitlines()
    out_lines = capsys.readouterr().out.splitlines()
    _check_logging_output(exp_lines, out_lines)


def test_viz_imshow_title_and_predefined_ax():
    dpm = DotPlotMatrix("AATCGATC", "TATCGAT", 6)
    t = "Hi I'm a title!"
    fig, ax = pyplot.subplots()
    viz_imshow(dpm, title=t, ax=ax)
    assert ax.get_title() == t
    # while we're at it, check that the x and y labels were correctly set
    assert ax.get_xlabel() == "$s_1$ (8 nt) \u2192"
    assert ax.get_ylabel() == "$s_2$ (7 nt) \u2192"


def test_viz_spy_title_and_predefined_ax():
    # like the above test, but with viz_spy()
    # sorry for the duplicated code, this is a bit sloppy
    dpm = DotPlotMatrix("AATCGATC", "TATCGAT", 6)
    t = "Hi I'm a title!"
    fig, ax = pyplot.subplots()
    viz_spy(dpm, title=t, ax=ax)
    assert ax.get_title() == t
    assert ax.get_xlabel() == "$s_1$ (8 nt) \u2192"
    assert ax.get_ylabel() == "$s_2$ (7 nt) \u2192"
