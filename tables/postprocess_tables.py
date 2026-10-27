"""Apply the manuscript's table terminology and rich-text notation rules.

This is the final presentation step after ``build_tables.py``. It changes only
listed text cells and their character-level typography. Numeric values, Boolean
values, formulas, original identifiers, cell styles and workbook structure are
checked against the input. S5 has no applicable changes and is copied unchanged.

Standalone usage (paths are independent of the current working directory)::

    python tables/postprocess_tables.py --input-dir tables/output
    python tables/postprocess_tables.py --input-dir input --output-dir output

The default operation is in place. ``formatting_changes.json`` records each
affected cell, the applied character ranges and preservation checks.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import re
import shutil
import tempfile
from collections import defaultdict
from zipfile import ZIP_DEFLATED, ZipFile

from lxml import etree

NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
X = "{" + NS["s"] + "}"
RID = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
XML_SPACE = "{http://www.w3.org/XML/1998/namespace}space"


def strings(archive):
    if "xl/sharedStrings.xml" not in archive.namelist():
        return []
    return [
        "".join(node.itertext())
        for node in etree.fromstring(archive.read("xl/sharedStrings.xml"))
    ]


def cell_text(cell, shared_strings):
    kind = cell.get("t")
    value = cell.find("s:v", NS)
    if kind == "inlineStr":
        return "".join(cell.find("s:is", NS).itertext())
    if kind == "s":
        return shared_strings[int(value.text)]
    if kind == "str":
        return value.text if value is not None else ""
    return None


def revise(index, sheet, address, text):
    """Reword only human-facing notes and column headers."""
    human = sheet == "说明" or re.fullmatch(r"[A-Z]+1", address)
    if not human:
        return text
    revised = re.sub(r"结瓜期叶(?!片)", "结瓜期叶片", text)
    revised = re.sub(r"苗期叶(?!片)", "苗期叶片", revised).replace("情境", "样品类型")
    revised = revised.replace("Context", "Sample type").replace("context", "sample type")
    revised = re.sub(r"(?<!相)同样品类型", "相同样品类型", revised)
    revised = revised.replace("该样品类型文库数", "该类样品文库数")
    revised = revised.replace("Libraries in sample type", "Libraries of this sample type")
    if index == 8:
        revised = revised.replace("一致性", "相关性").replace("concordance", "correlation")
        revised = re.sub(r"Included in correlation(?! analysis)", "Included in correlation analysis", revised)
    if index == 3 and sheet == "说明":
        revised = revised.replace("达到门槛", "达到阈值").replace("全基因水平的DESeq2结果", "全部受检基因的DESeq2结果")
    if index == 4 and sheet == "说明":
        revised = revised.replace("果实取样时期间", "果实取样时期之间").replace("两个取样时期间", "两个取样时期之间")
        revised = revised.replace("正值表示偏离随取样时期增大", "正值表示较晚取样时期的表达偏离更大")
        revised = revised.replace("between fruit stages", "between fruit sampling times")
    if index == 7 and sheet == "说明":
        revised = revised.replace("支持门槛", "支持文库数阈值")
        revised = revised.replace("D为相同样品类型", "D为同一样品类型中")
    # These underscores are mathematical labels in prose, never raw identifiers.
    if index == 2 and sheet == "说明" and address == "A14":
        revised = revised.replace("W_MP", "WMP").replace("N_MP", "NMP")
        revised = revised.replace("c_F1", "cF₁").replace("c_C", "cC").replace("c_P", "cP")
    if "D_MD" in revised or "D_MS" in revised:
        revised = revised.replace("D_MD", "DMD").replace("D_MS", "DMS")
    return revised


def spans(index, sheet, address, text):
    """Return ordered typography rules; later rules override overlapping ones."""
    output = []
    human = sheet == "说明" or bool(re.fullmatch(r"[A-Z]+1", address))
    categorical = index == 8 and sheet == "跨平台效应匹配" and text == "F1–parent log2FC"
    gene_label = index == 3 and sheet == "图3C基因" and re.fullmatch(r"B(?:[2-9]|1[0-9]|2[0-4])", address)
    reference_gene = index == 8 and (
        sheet == "引物序列" and text == "Actin2" or sheet == "说明" and "Actin2" in text
    )
    if not (human or categorical or gene_label or reference_gene):
        return output

    def add(start, end, **properties):
        output.append({"start": start, "end": end, "text": text[start:end], **properties})

    for match in re.finditer(r"(?<![A-Za-z0-9_])F1(?![A-Za-z0-9_])", text):
        add(match.start() + 1, match.end(), vertAlign="subscript", italic=False)
    for match in re.finditer(r"log(2|10)", text):
        add(match.start() + 3, match.end(), vertAlign="subscript", italic=False)
    for match in re.finditer(r"(?<![A-Za-z])Ct(?![A-Za-z])", text):
        add(match.start() + 1, match.end(), vertAlign="subscript", italic=False)
        add(match.start(), match.start() + 1, italic=True)
    for match in re.finditer(r"df([12])", text):
        add(match.start() + 2, match.end(), vertAlign="subscript", italic=False)
    if reference_gene:
        for match in re.finditer("Actin2", text):
            add(match.start(), match.end(), italic=True)
    if gene_label:
        match = re.match(r"[A-Z][A-Z0-9]*", text)
        if match:
            add(match.start(), match.end(), italic=True)
    if human:
        for match in re.finditer(r"(?<![A-Za-z0-9_])[Pq](?![A-Za-z0-9_])", text):
            if match.group() == "P" and (
                text[match.end():].startswith("在MD")
                or match.start() > 0 and text[match.start() - 1] == "μ"
            ):
                continue
            add(match.start(), match.end(), italic=True)
        for match in re.finditer(r"(?<![A-Za-z0-9_])h(?![A-Za-z0-9_])|(?<![A-Za-z])t(?=检验)", text):
            add(match.start(), match.end(), italic=True)
        for match in re.finditer(r"ρ|μ", text):
            add(match.start(), match.end(), italic=True)
        for match in re.finditer(r"ΔR|ΔD", text):
            add(match.end() - 1, match.end(), italic=True)
        for match in re.finditer(r"(?<=偏离)D|(?<=\|)D(?=\|)|\bD(?=\s*[=(])|^D(?=为)|(?<=Departure )D|(?<![A-Za-z0-9_])D(?=的q|$)|(?<=[（(])D(?=[）)])", text):
            add(match.start(), match.end(), italic=True)
        for match in re.finditer(r"(?<![A-Za-z0-9_])F(?=值|\s|$)", text):
            add(match.start(), match.end(), italic=True)
        for match in re.finditer(r"D(MD|MS)", text):
            add(match.start(), match.start() + 1, italic=True)
            add(match.start() + 1, match.end(), vertAlign="subscript", italic=False)
        for match in re.finditer(r"μ(F1|F₁|C|P)", text):
            add(match.start() + 1, match.end(), vertAlign="subscript", italic=False)
        if index == 7 and sheet == "说明" and address == "A7":
            for match in re.finditer(r"(?<=\))(MS|MD)", text):
                add(match.start(), match.end(), vertAlign="subscript", italic=False)
        if index == 2 and sheet == "说明" and address == "A14":
            for match in re.finditer(r"(?:N̄|W̄)(C|P)|[WN](MP)|δ([NW])|c(F₁|F1|C|P)", text):
                for group in range(1, 5):
                    if match.group(group):
                        add(match.start(group), match.end(group), vertAlign="subscript", italic=match.group(group) in ["N", "W"])
                        add(match.start(), match.start(group), italic=True)
                        break
    return output


def sheet_parts(archive):
    workbook = etree.fromstring(archive.read("xl/workbook.xml"))
    relationships = {
        node.get("Id"): node.get("Target")
        for node in etree.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    }
    for sheet in workbook.findall("s:sheets/s:sheet", NS):
        target = relationships[sheet.get(RID)]
        yield sheet.get("name"), target.lstrip("/") if target.startswith("/") else "xl/" + target


def all_cells(archive):
    shared_strings = strings(archive)
    output = {}
    for _, part in sheet_parts(archive):
        root = etree.fromstring(archive.read(part))
        for cell in root.findall(".//s:c", NS):
            output[(part, cell.get("r"))] = {
                "type": cell.get("t"),
                "value": cell.findtext("s:v", namespaces=NS),
                "formula": cell.findtext("s:f", namespaces=NS),
                "text": cell_text(cell, shared_strings),
                "style": cell.get("s", "0"),
            }
    return output


def base_rpr(font):
    properties = etree.Element(X + "rPr")
    for child in font:
        node = copy.deepcopy(child)
        if etree.QName(node).localname == "name":
            node.tag = X + "rFont"
        properties.append(node)
    return properties


def richtext(cell, text, ranges, font):
    properties = [{} for _ in text]
    for span in ranges:
        assert text[span["start"]:span["end"]] == span["text"]
        for index in range(span["start"], span["end"]):
            properties[index].update({k: span[k] for k in ("italic", "vertAlign") if k in span})
    for child in list(cell):
        if etree.QName(child).localname in ("v", "is"):
            cell.remove(child)
    cell.set("t", "inlineStr")
    inline = etree.SubElement(cell, X + "is")
    if not ranges:
        node = etree.SubElement(inline, X + "t")
        node.text = text
        if text.strip() != text:
            node.set(XML_SPACE, "preserve")
        return
    boundaries = [0] + [i for i in range(1, len(text)) if properties[i] != properties[i - 1]] + [len(text)]
    order = {name: i for i, name in enumerate(("rFont", "charset", "family", "b", "i", "strike", "outline", "shadow", "condense", "extend", "color", "sz", "u", "vertAlign", "scheme"))}
    for start, end in zip(boundaries, boundaries[1:]):
        run = etree.SubElement(inline, X + "r")
        run_properties = base_rpr(font)
        for key, value in properties[start].items():
            tag = "i" if key == "italic" else key
            for old in run_properties.findall("s:" + tag, NS):
                run_properties.remove(old)
            node = etree.SubElement(run_properties, X + tag)
            node.set("val", ("1" if value else "0") if key == "italic" else value)
        for child in sorted(list(run_properties), key=lambda item: order.get(etree.QName(item).localname, 99)):
            run_properties.remove(child)
            run_properties.append(child)
        run.append(run_properties)
        node = etree.SubElement(run, X + "t")
        node.text = text[start:end]
        if node.text.strip() != node.text:
            node.set(XML_SPACE, "preserve")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def process_workbook(source, output, index):
    source_hash = sha256(source)
    edits = []
    with ZipFile(source) as archive:
        before = all_cells(archive)
        shared_strings = strings(archive)
        for sheet, part in sheet_parts(archive):
            root = etree.fromstring(archive.read(part))
            for cell in root.findall(".//s:c", NS):
                address = cell.get("r")
                text = cell_text(cell, shared_strings)
                if text is None:
                    continue
                revised = revise(index, sheet, address, text)
                ranges = spans(index, sheet, address, revised)
                if revised != text or ranges:
                    edits.append({"sheet": sheet, "part": part, "cell": address, "before": text, "after": revised, "spans": ranges, "is_value_change": revised != text})
        updates = {}
        if edits:
            styles = etree.fromstring(archive.read("xl/styles.xml"))
            fonts = styles.find("s:fonts", NS)
            cell_styles = styles.find("s:cellXfs", NS)
            grouped = defaultdict(list)
            for edit in edits:
                grouped[edit["part"]].append(edit)
            for part, changes in grouped.items():
                root = etree.fromstring(archive.read(part))
                original_root = copy.deepcopy(root)
                for edit in changes:
                    cell = root.find(f'.//s:c[@r="{edit["cell"]}"]', NS)
                    font_index = int(cell_styles[int(cell.get("s", "0"))].get("fontId", "0"))
                    richtext(cell, edit["after"], edit["spans"], fonts[font_index])
                # Restoring the listed cells must restore the whole sheet XML.
                restored = copy.deepcopy(root)
                for edit in changes:
                    current = restored.find(f'.//s:c[@r="{edit["cell"]}"]', NS)
                    original = original_root.find(f'.//s:c[@r="{edit["cell"]}"]', NS)
                    current.getparent().replace(current, copy.deepcopy(original))
                assert etree.tostring(restored, method="c14n") == etree.tostring(original_root, method="c14n")
                encoded = etree.tostring(root, encoding="utf-8")
                if encoded != archive.read(part):
                    updates[part] = encoded
        output.parent.mkdir(parents=True, exist_ok=True)
        if updates:
            with tempfile.NamedTemporaryFile(dir=output.parent, suffix=".xlsx", delete=False) as handle:
                temporary = Path(handle.name)
            try:
                with ZipFile(temporary, "w", ZIP_DEFLATED) as result:
                    for info in archive.infolist():
                        result.writestr(copy.copy(info), updates.get(info.filename, archive.read(info.filename)))
                with ZipFile(temporary) as result:
                    after = all_cells(result)
                    assert before.keys() == after.keys()
                    expected = {(e["part"], e["cell"]): e["after"] for e in edits}
                    for key, old in before.items():
                        new = after[key]
                        assert old["style"] == new["style"] and old["formula"] == new["formula"], key
                        if old["text"] is None:
                            assert old == new, key
                        else:
                            assert new["text"] == expected.get(key, old["text"]), key
                    changed = [name for name in archive.namelist() if archive.read(name) != result.read(name)]
                    assert set(changed) == set(updates)
                temporary.replace(output)
            finally:
                if temporary.exists():
                    temporary.unlink()
        elif source.resolve() != output.resolve():
            shutil.copy2(source, output)
    nontext = [cell for cell in before.values() if cell["text"] is None and cell["value"] is not None]
    result = {
        "file": source.name,
        "text_cells_changed": sum(e["is_value_change"] for e in edits),
        "format_cells": sum(bool(e["spans"]) for e in edits),
        "nontext_value_cells_preserved": len(nontext),
        "numeric_cells_preserved": sum(cell["type"] in (None, "n") for cell in nontext),
        "boolean_cells_preserved": sum(cell["type"] == "b" for cell in nontext),
        "formula_cells_preserved": sum(cell["formula"] is not None for cell in before.values()),
        "changed_parts": sorted(updates),
        "other_zip_parts_byte_identical": True,
        "noncell_xml_unchanged": True,
        "input_sha256": source_hash,
        "output_sha256": sha256(output),
        "edits": edits,
    }
    if index == 5:
        assert result["input_sha256"] == result["output_sha256"], "S5 must remain byte-identical"
    return result


def postprocess_tables(input_dir, output_dir=None, manifest=None):
    input_dir = Path(input_dir).resolve()
    output_dir = Path(output_dir).resolve() if output_dir is not None else input_dir
    inputs = [input_dir / f"S{index}.xlsx" for index in range(1, 9)]
    missing = [path.name for path in inputs if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing table workbooks: " + ", ".join(missing))
    results = [process_workbook(path, output_dir / path.name, index) for index, path in enumerate(inputs, 1)]
    manifest = Path(manifest) if manifest is not None else output_dir / "formatting_changes.json"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"text_cells_changed": sum(r["text_cells_changed"] for r in results), "format_cells": sum(r["format_cells"] for r in results), "S5_byte_identical": results[4]["input_sha256"] == results[4]["output_sha256"]}, ensure_ascii=False))
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path(__file__).resolve().parent / "output")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--manifest", type=Path)
    arguments = parser.parse_args()
    postprocess_tables(arguments.input_dir, arguments.output_dir, arguments.manifest)
