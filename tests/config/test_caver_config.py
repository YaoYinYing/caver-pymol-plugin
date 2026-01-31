import json
from textwrap import dedent

import pytest

from Caver4 import caver_config
from Caver4.caver_config import CaverConfig


def test_all_keys_hide_ui_and_start_point_fields():
    config = CaverConfig()
    hidden = {"output_dir", "selection_name", "start_point_x", "start_point_y", "start_point_z"}
    result = set(config.all_keys)
    assert hidden.isdisjoint(result)
    assert {"shell_radius", "probe_radius", "number_of_approximating_balls"} <= result


def test_has_get_set_and_delete_cycle():
    config = CaverConfig()
    config.set("probe_radius", 2.3)
    assert config.has("probe_radius")
    assert config.get("probe_radius") == 2.3
    config.delete("probe_radius")
    assert not config.has("probe_radius")


def test_from_json_casts_types_and_preserves_extra_keys(tmp_path):
    json_payload = {
        "shell_depth": "6",
        "shell_radius": "4.2",
        "number_of_approximating_balls": 7,
        "custom_path": ["a", "b"],
    }
    json_file = tmp_path / "config.json"
    json_file.write_text(json.dumps(json_payload))

    config = CaverConfig.from_json(str(json_file))

    assert config.shell_depth == 6
    assert pytest.approx(config.shell_radius) == 4.2
    assert config.number_of_approximating_balls == 7
    assert config.custom_path == ["a", "b"]


def test_to_json_persists_entire_state(tmp_path):
    config = CaverConfig(shell_depth=9)
    config.extra_flag = "enabled"
    target = tmp_path / "out.json"

    config.to_json(str(target))
    written = json.loads(target.read_text())

    assert written["shell_depth"] == 9
    assert written["extra_flag"] == "enabled"
    assert written["probe_radius"] == pytest.approx(config.probe_radius)


def test_from_txt_interprets_comments_and_values(tmp_path):
    txt_payload = dedent(
        """
        # comment
        shell_depth 10 # trailing
        shell_radius 5.5
        custom_text "value with spaces"
        orphan_key
        """
    ).strip()
    text_file = tmp_path / "config.txt"
    text_file.write_text(txt_payload)

    config = CaverConfig.from_txt(str(text_file))

    assert config.shell_depth == 10
    assert config.shell_radius == pytest.approx(5.5)
    assert config.custom_text == '"value with spaces"'
    assert not hasattr(config, "orphan_key")


@pytest.mark.parametrize(
    "attr_name,new_value,expected",
    [
        ("shell_radius", "5.1", 5.1),
        ("shell_depth", "8", 8),
        ("number_of_approximating_balls", "12", 12),
    ],
)
def test__set_value_casts_numeric_types(attr_name, new_value, expected):
    config = CaverConfig()
    config._set_value(attr_name, new_value)
    assert getattr(config, attr_name) == expected


def test__set_value_handles_boolean_and_unknown_strings():
    config = CaverConfig()
    config.custom_toggle = False
    config._set_value("custom_toggle", "yes")
    config._set_value("new_field", "raw")
    assert config.custom_toggle is True
    assert config.new_field == "raw"


def test__get_value_for_boolean_and_float():
    config = CaverConfig()
    config.boolean_flag = True
    config.float_flag = 1.25
    assert config._get_value("boolean_flag") == "yes"
    assert config._get_value("float_flag") == "1.25"


@pytest.mark.parametrize("value,expected", [(True, "yes"), (False, "no")])
def test__yes_or_no(value, expected):
    assert CaverConfig._yes_or_no(value) == expected


@pytest.mark.parametrize("text,expected", [("yes", True), ("no", False), ("anything", False)])
def test__true_or_false(text, expected):
    assert CaverConfig._true_or_false(text) is expected


def test_to_txt_uses_template_and_formats_values(tmp_path, monkeypatch):
    template = dedent(
        """
        # used as a template
        # header

        shell_radius 0
        custom_toggle yes
        float_value 0.0
        starting_point_coordinates ???
        str_field ???
        """
    ).strip()
    template_file = tmp_path / "template.txt"
    template_file.write_text(template)
    monkeypatch.setattr(caver_config, "CONFIG_TXT", str(template_file))

    config = CaverConfig(shell_radius=4.567)
    config.custom_toggle = True
    config.float_value = 1.23
    config.str_field = "???"
    config.start_point_x = 1
    config.start_point_y = 2
    config.start_point_z = 3

    out_path = tmp_path / "export" / "config.txt"
    config.to_txt(str(out_path))
    lines = out_path.read_text().splitlines()

    assert lines == [
        "# header",
        "",
        "shell_radius 4.6",
        "custom_toggle yes",
        "float_value 1.2",
        "starting_point_coordinates 1 2 3",
    ]
