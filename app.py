import json
import base64
import posixpath
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape as xml_escape
import re
import hashlib
import zipfile
from pathlib import Path
from io import StringIO, BytesIO

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import font_manager

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from Bio import SeqIO, Restriction

MM = 1.0 / 25.4  # mm → inch

st.set_page_config(layout="wide", page_title="논문용 유전자 맵 스튜디오")
st.title("🧬 논문 피규어용 유전자 맵 스튜디오")

with st.expander("🔒 데이터 보안 및 면책 조항 (사용 전 필독)", expanded=False):
    st.markdown("""
**📁 데이터 처리 방식**
본 프로그램은 사용자의 데이터를 수집하거나 외부 서버에 저장하지 않습니다.
업로드된 GenBank 파일은 브라우저 세션이 유지되는 동안 **서버 메모리에서만** 처리되며,
페이지 새로고침·탭 종료 시 즉시 파기됩니다. 미발표 서열, 특허 출원 전 데이터,
환자 유래 민감 정보 등도 안심하고 사용하실 수 있습니다.

**⚠️ 면책 조항 (Disclaimer)**
본 소프트웨어는 **'있는 그대로(As-is)'** 제공되며, 어떠한 형태의 명시적·묵시적 보증도 하지 않습니다.
프로그램은 좌표 변환, 자동 분류 등의 과정에서 예기치 않은 오류를 포함할 수 있으며,
**생성된 피규어의 좌표·라벨·결과물에 대한 최종 검증 책임은 전적으로 사용자 본인에게 있습니다.**
논문 투고 전 반드시 원본 GenBank 파일과 대조하여 정확성을 확인해 주세요.

**📜 License & Author / 라이선스 및 제작**
This software was developed by **YJ LEE (Yujin Lee)** and is released
for **NONCOMMERCIAL ACADEMIC USE ONLY** (PolyForm Noncommercial 1.0.0
based). Academic, educational, and nonprofit research use is freely
permitted. Any commercial use is strictly prohibited without prior
written permission from the author.

본 소프트웨어는 **YJ LEE (Yujin Lee)** 가 제작하였으며,
**비영리 학술 사용 전용** (PolyForm Noncommercial 1.0.0 기반)으로
배포됩니다. 학술·교육·비영리 연구 목적의 사용은 자유롭게 허용되나,
상업적 사용은 저작자의 사전 서면 허가 없이 엄격히 금지됩니다.
상업적 사용을 원하시는 경우 저작자에게 별도 문의 바랍니다.
    """)


# ============================================================
# 폰트 설정 — Arial 우선, 저널 투고용 벡터 옵션 포함
# ============================================================
@st.cache_resource
def setup_font():
    available = {f.name for f in font_manager.fontManager.ttflist}
    chosen = "DejaVu Sans"
    for candidate in ["Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "Arimo"]:
        if candidate in available:
            chosen = candidate
            break

    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = [chosen, "Liberation Sans", "DejaVu Sans"]
    plt.rcParams["mathtext.default"] = "regular"
    try:
        plt.rcParams["mathtext.fontset"] = "custom"
        plt.rcParams["mathtext.rm"] = chosen
        plt.rcParams["mathtext.it"] = f"{chosen}:italic"
        plt.rcParams["mathtext.bf"] = f"{chosen}:bold"
        plt.rcParams["mathtext.bfit"] = f"{chosen}:italic:bold"
    except Exception:
        plt.rcParams["mathtext.fontset"] = "stixsans"

    # 투고용: PDF/EPS에 폰트를 TrueType으로 embed, SVG는 텍스트를 살려서 내보냄
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"
    return chosen


ACTIVE_FONT = setup_font()


# ============================================================
# 프리셋 기본값
# ============================================================
PRESET_DEFAULTS = {
    # 텍스트
    "wt_label": "WT allele",
    "mut_label": "Mutant allele",
    "text_color": "#000000",
    "label_font_size": 11,
    "marker_font_size": 9,
    "enzyme_font_size": 9,
    "primer_font_size": 9,
    # 색상
    "color_exon": "#002060", "color_marker": "#FFC000",
    "color_enzyme": "#FF0000", "color_promoter": "#AAAAAA",
    "color_probe": "#92D050", "color_line": "#000000",
    # 캔버스 (가로만 지정, 세로는 내용에 맞춰 자동)
    "fig_width": 7.2,
    "label_area_mm": 26.0,
    "dpi_setting": 600,
    # 요소 크기 (mm)
    "box_h_mm": 3.0,
    "probe_h_mm": 3.0,
    "enzyme_tick_mm": 2.5,
    "primer_stem_mm": 3.5,
    "arrow_mm": 2.5,
    "track_gap_mm": 6.0,
    "line_width": 1.5,
    # 범위
    "trim_left": 0, "trim_right": 0,
    # 옵션
    "show_scalebar": True, "scalebar_bp": 1000,
    "scalebar_color": "#000000", "scalebar_lw": 2.5,
    "show_hr": True,
    "dedup_enzyme": False,
    "exon_arrow": True,
    "exon_head_ratio": 1.8,
    "auto_stagger": True,
    # 자동 인식 규칙 (직접 수정 가능)
    "rule_marker": "NAT, NEO, HYG, HPH, TRP1, URA5, URA3, ADE2, G418, BLE, SAT1, KAN, ZEO, NOU, PUR",
    "rule_promoter": "PROMOTER, PROM",
    "rule_probe": "PROBE",
    "rule_primer": "PRIMER, FWD, REV",
    "rule_enzyme_max": 12,
    "rule_primer_max": 60,
    "hide_unnamed": True,
}

CATEGORIES = ["Exon", "Marker", "Promoter", "Primer", "Enzyme", "Probe"]
ARROW_CHOICES = ["자동", "없음", "→", "←"]

for _k, _v in PRESET_DEFAULTS.items():
    if _k not in st.session_state:
        st.session_state[_k] = _v


# ============================================================
# 사용 가이드 Word 생성 — 별도 DOCX 파일이나 추가 패키지 없이 다운로드
# ============================================================
@st.cache_data
def generate_usage_guide():
    def paragraph(text, style="Normal"):
        return ('<w:p><w:pPr><w:pStyle w:val="' + style + '"/></w:pPr>'
                '<w:r><w:t xml:space="preserve">' + xml_escape(text) + '</w:t></w:r></w:p>')
    def heading(text):
        result = paragraph(text, "Heading1")
        if text in {"엑셀 입력 양식", "여러 디자인을 한 번에 출력하기"}:
            result = result.replace("<w:pPr>", "<w:pPr><w:pageBreakBefore/>", 1)
        return result
    def table(headers, rows, widths):
        cells = [headers] + rows
        parts = ['<w:tbl><w:tblPr><w:tblW w:w="9972" w:type="dxa"/><w:tblBorders>']
        for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
            parts.append('<w:' + edge + ' w:val="single" w:sz="4" w:color="D9D9D9"/>')
        parts.append('</w:tblBorders><w:tblCellMar><w:top w:w="90" w:type="dxa"/><w:left w:w="120" w:type="dxa"/><w:bottom w:w="90" w:type="dxa"/><w:right w:w="120" w:type="dxa"/></w:tblCellMar></w:tblPr><w:tblGrid>')
        for width in widths: parts.append('<w:gridCol w:w="' + str(width) + '"/>')
        parts.append('</w:tblGrid>')
        for index, row in enumerate(cells):
            parts.append('<w:tr>' + ('<w:trPr><w:tblHeader/></w:trPr>' if index == 0 else ''))
            for width, text in zip(widths, row):
                fill = '233D59' if index == 0 else ('EDF3FA' if index % 2 == 1 else 'FFFFFF')
                parts.append('<w:tc><w:tcPr><w:tcW w:w="' + str(width) + '" w:type="dxa"/><w:shd w:fill="' + fill + '"/><w:vAlign w:val="center"/></w:tcPr><w:p><w:pPr><w:spacing w:after="40"/></w:pPr><w:r><w:rPr><w:sz w:val="20"/>' + ('<w:b/><w:color w:val="FFFFFF"/>' if index == 0 else '') + '</w:rPr><w:t>' + xml_escape(str(text)) + '</w:t></w:r></w:p></w:tc>')
            parts.append('</w:tr>')
        parts.append('</w:tbl>')
        return ''.join(parts)
    page_break = ''
    body = [paragraph("유전자 맵 스튜디오 사용 가이드", "Title"),
        paragraph("WT와 Mutant의 GenBank 주석을 나란히 그려 PNG, SVG, PDF로 저장하는 도구입니다. 한 디자인을 직접 편집하거나, 엑셀에 파일 짝과 라벨 및 단어별 기울임을 입력하여 여러 그림을 일괄 출력할 수 있습니다."),
        heading("한 디자인 만들기"),
        paragraph("1. 파일 업로드 탭에서 WT와 Mutant의 .gb 또는 .gbk 파일을 각각 선택합니다."),
        paragraph("2. 요소 편집 탭에서 표시 여부, 이름, 종류, 화살표 방향, 위아래 위치를 확인합니다. 화면에서 표시하는 좌표는 업로드 파일의 주석을 따릅니다."),
        paragraph("3. HR 연결 탭에서 점선으로 연결할 이름을 지정합니다. 기본값은 L1과 L2를 서로 교차하여 연결하고 R1과 R2는 같은 이름끼리 연결합니다. 파일에 없는 이름은 연결하지 않습니다."),
        paragraph("4. 미리보기 및 다운로드 탭에서 PNG, SVG, PDF 중 원하는 형식으로 저장합니다."),
        heading("라벨과 글씨 서식"),
        paragraph("왼쪽 이름의 기본값은 WT allele과 Mutant allele입니다. 사이드바에서 직접 바꾸거나 엑셀에서 파일 짝마다 지정할 수 있습니다. 엑셀에 일치하는 파일 짝이 있으면 엑셀의 라벨을 사용합니다."),
        paragraph("수동 입력에서는 *기울임*, **굵게**, ***굵은기울임***을 사용합니다. ~아래첨자~와 ^위첨자^도 지원합니다. 예를 들어 P~CTR4~: MSS4는 CTR4를 아래첨자로 표시합니다."),
        heading("박스와 이름 배치"),
        paragraph("Marker는 각 주석의 strand 방향으로 화살표가 표시됩니다. 요소 편집에서 없음, 오른쪽 또는 왼쪽으로 바꿀 수 있습니다. 마지막 exon도 strand 방향의 화살표로 표시합니다."),
        paragraph("Mutant에서 박스 안에 들어가지 않는 긴 Marker와 Promoter 이름은 박스 바로 아래에 표시합니다. 모든 Probe 이름은 박스의 폭과 관계없이 바로 아래에 표시합니다. 제한효소 눈금은 시퀀스 선에 붙여 둡니다."),
        heading("크기와 색상"),
        paragraph("사이드바에서 그림 너비, 왼쪽 라벨 공간, 글씨 크기, 색상, 박스 높이, 트랙 간격과 scale bar를 조절합니다. PNG의 DPI는 래스터 해상도입니다. SVG와 PDF는 벡터 형식으로 저장됩니다."),
        page_break,
        heading("엑셀 입력 양식"),
        paragraph("엑셀 설정 탭에서 입력 양식을 내려받습니다. 예시 행의 파일명을 실제 파일명으로 바꾸고 Designs와 TextStyles를 작성한 뒤 같은 탭에 업로드합니다. 시트 이름과 첫 행의 열 이름은 그대로 유지합니다."),
        heading("Designs 시트"),
        table(["열", "입력 내용"], [
            ["디자인 ID", "D001 등의 고유 ID. TextStyles와 연결하는 기준입니다."],
            ["출력", "예는 출력, 아니오는 일괄 작업에서 제외합니다."],
            ["WT 파일 / Mutant 파일", "GenBank 파일명을 확장자까지 정확하게 입력합니다."],
            ["출력 이름", "ZIP 안의 폴더와 그림 파일 이름을 지정합니다."],
            ["WT 라벨 / Mutant 라벨", "각 시퀀스 왼쪽에 표시할 이름을 지정합니다. 빈칸은 기본 라벨로 채웁니다."]], [2900,7072]),
        heading("TextStyles 시트"),
        table(["열", "입력 내용"], [
            ["디자인 ID", "*는 모든 디자인, D001 등의 ID는 해당 디자인만 적용합니다."],
            ["적용 대상", "전체, WT 라벨, Mutant 라벨, WT 내부, Mutant 내부 중 선택합니다."],
            ["단어", "서식을 바꿀 글씨를 한 행에 하나씩 입력합니다."],
            ["기울임", "예는 기울임, 아니오는 일반체입니다. 기존의 기울임도 해제할 수 있습니다."],
            ["일치 방식", "부분 일치, 단어 일치, 전체 일치 중 선택합니다."]], [2900,7072]),
        paragraph("부분 일치는 글씨 안에 포함된 부분만 바꿉니다. 단어 일치는 독립된 단어에 적용하며, 전체 일치는 글씨 전체가 같은 경우에만 적용합니다. 대소문자는 구분하지 않습니다."),
        paragraph("예: 디자인 ID D001, 적용 대상 전체, 단어 CTR4, 기울임 예, 일치 방식 부분 일치로 입력하면 CTR4 promoter에서 CTR4만 기울어집니다. 아래첨자 CTR4에도 적용됩니다."),
        paragraph("여러 규칙이 같은 글씨에 겹치면 TextStyles의 아래쪽 행이 우선합니다. 기울임 선택만 바뀌며 기존 굵기와 위아래첨자는 유지됩니다. $...$로 직접 입력한 수식에는 단어별 규칙을 적용하지 않습니다."),
        page_break,
        heading("여러 디자인을 한 번에 출력하기"),
        paragraph("1. 엑셀 설정 탭에 작성한 양식을 업로드합니다. 양식 없이도 일괄 작업을 시작할 수 있습니다."),
        paragraph("2. 일괄 작업 탭에서 WT 파일과 Mutant 파일을 각각 여러 개 선택합니다."),
        paragraph("3. 표에서 파일 짝과 라벨을 확인합니다. 엑셀을 업로드했다면 그 행을 사용합니다. 엑셀이 없으면 파일명의 CNAG 번호 또는 공통 이름으로 짝을 추천하며, 모호한 짝은 직접 지정합니다. 같은 WT를 여러 Mutant와 짝지을 수 있습니다."),
        paragraph("4. 필요한 출력 행과 PNG, SVG, PDF 형식을 선택합니다. 사이드바의 그림 설정과 일괄 HR 연결이 모든 디자인에 적용됩니다."),
        paragraph("5. 선택한 디자인 일괄 생성을 누른 뒤 ZIP을 다운로드합니다. 디자인마다 폴더가 생성되고 처리 결과는 processing_report.csv에 기록됩니다. 일부 파일에서 오류가 나도 나머지 파일은 처리합니다."),
        heading("한 파일 짝에 엑셀 적용하기"),
        paragraph("파일 업로드 탭의 WT와 Mutant가 엑셀의 파일명과 일치하면 해당 라벨과 TextStyles가 자동으로 적용됩니다. 같은 파일 짝을 여러 디자인에 등록했다면 적용할 디자인 ID를 선택합니다."),
        heading("직접 편집과 저장"),
        paragraph("직접 편집 탭에서 글씨, 눈금과 박스를 클릭하여 드래그하거나 방향키로 이동할 수 있습니다. Shift를 함께 누르면 더 크게 이동합니다. 설정을 바꾸면 그림이 다시 생성되므로 세부 위치 조정은 마지막에 합니다."),
        paragraph("직접 편집한 그림은 그 탭의 SVG 저장 또는 PNG 저장 버튼으로 내려받습니다. 직접 편집 내용은 미리보기 탭의 원본 PDF나 일괄 출력 ZIP에 반영되지 않습니다."),
        heading("글꼴과 문제 해결"),
        paragraph("서버에 Arial이 없으면 대체 글꼴을 사용합니다. Streamlit Community Cloud에서는 GitHub 저장소 최상위에 packages.txt를 만들고 fonts-liberation을 한 줄로 입력하면 Liberation Sans를 사용할 수 있습니다."),
        paragraph("파일 짝이 적용되지 않으면 엑셀과 업로드 파일의 이름 및 확장자를 대조합니다. 양식에는 수식 대신 텍스트를 입력합니다. 긴 이름은 줄이거나 그림 너비와 왼쪽 라벨 공간을 늘려 표시를 조정합니다."),
        heading("사용 범위"),
        paragraph("제작 YJ LEE Yujin Lee. 비영리 학술 사용 전용입니다. 최종 그림은 업로드한 원본 GenBank의 주석과 대조하여 확인합니다. 이 가이드는 프로그램에서 생성되므로 별도의 Word 파일을 배포 저장소에 올릴 필요가 없습니다.")]
    styles = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:docDefaults><w:rPrDefault><w:rPr><w:rFonts w:ascii="Arial" w:hAnsi="Arial" w:eastAsia="맑은 고딕"/><w:sz w:val="22"/><w:color w:val="000000"/></w:rPr></w:rPrDefault><w:pPrDefault><w:pPr><w:spacing w:after="100" w:line="252" w:lineRule="auto"/></w:pPr></w:pPrDefault></w:docDefaults><w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/></w:style><w:style w:type="paragraph" w:styleId="Title"><w:name w:val="Title"/><w:basedOn w:val="Normal"/><w:pPr><w:keepNext/><w:spacing w:after="220"/></w:pPr><w:rPr><w:b/><w:sz w:val="36"/><w:color w:val="000000"/></w:rPr></w:style><w:style w:type="paragraph" w:styleId="Heading1"><w:name w:val="heading 1"/><w:basedOn w:val="Normal"/><w:pPr><w:keepNext/><w:spacing w:before="140" w:after="80"/><w:outlineLvl w:val="0"/></w:pPr><w:rPr><w:b/><w:sz w:val="26"/><w:color w:val="000000"/></w:rPr></w:style></w:styles>'''
    document = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:body>' + ''.join(body) + '<w:sectPr><w:pgSz w:w="12240" w:h="15840"/><w:pgMar w:top="900" w:right="1134" w:bottom="900" w:left="1134" w:header="0" w:footer="0"/></w:sectPr></w:body></w:document>'
    content_types = '<?xml version="1.0" encoding="UTF-8"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/><Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/></Types>'
    root_rels = '<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/></Relationships>'
    word_rels = '<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>'
    buffer = BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, xml in [("[Content_Types].xml", content_types), ("_rels/.rels", root_rels),
                          ("word/document.xml", document), ("word/styles.xml", styles),
                          ("word/_rels/document.xml.rels", word_rels)]:
            archive.writestr(name, xml.encode("utf-8"))
    return buffer.getvalue()


# ============================================================
# 사이드바
# ============================================================
with st.sidebar:
    st.header("⚙️ 컨트롤 패널")
    st.caption(f"현재 사용 중인 글꼴: **{ACTIVE_FONT}**")
    if ACTIVE_FONT not in ("Arial", "Helvetica", "Liberation Sans", "Arimo"):
        st.warning(
            "서버에 Arial이 설치되어 있지 않아 대체 글꼴로 그려집니다. "
            "`packages.txt`에 `fonts-liberation`을 추가하면 Arial과 자간·자폭이 같은 "
            "Liberation Sans가 적용됩니다."
        )

    with st.expander("ℹ️ 사용 가이드 다운로드", expanded=False):
        st.caption("유전자 맵 스튜디오의 상세 사용법과 꿀팁을 확인하세요.")
        st.download_button(
            label="📄 가이드 다운로드 (.docx)", data=generate_usage_guide(),
            file_name="Gene_Map_Studio_User_Guide.docx",
            mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
            use_container_width=True,
        )

    with st.expander("📐 크기", expanded=True):
        st.caption(
            "가로 길이만 정하면 세로는 요소 개수에 맞춰 자동으로 계산됩니다. "
            "가로를 늘리면 글자 크기는 그대로인 채 유전자 영역만 넓게 펼쳐집니다. "
            "저널 기준: 1단 약 3.5 in, 2단 약 7.2 in."
        )
        fig_width = st.slider("가로 길이 (inch)", 4.0, 20.0, step=0.5, key="fig_width")
        label_area_mm = st.slider("왼쪽 라벨 공간 (mm)", 0.0, 60.0, step=2.0, key="label_area_mm")
        dpi_setting = st.select_slider("PNG 해상도 (DPI)", options=[150, 300, 600, 1200], key="dpi_setting")

    with st.expander("📝 라벨 & 텍스트", expanded=True):
        st.caption(
            "`**굵게**` · `*기울임*` · `^위첨자^` · `~아래첨자~`\n\n"
            "**조합 가능:** `***굵은기울임***`, `**car1^Δ^**`, `*p~TEF1~*`"
        )
        wt_label = st.text_input("WT 라벨", key="wt_label")
        mut_label = st.text_input("Mutant 라벨", key="mut_label")
        text_color = st.color_picker("텍스트 색상", key="text_color")
        label_font_size = st.slider("라벨 (WT/Mut)", 6, 24, key="label_font_size")
        marker_font_size = st.slider("Marker / Promoter", 6, 20, key="marker_font_size")
        enzyme_font_size = st.slider("Enzyme", 6, 20, key="enzyme_font_size")
        primer_font_size = st.slider("Primer", 6, 20, key="primer_font_size")

    with st.expander("🎨 색상"):
        color_exon = st.color_picker("Exon", key="color_exon")
        color_marker = st.color_picker("Marker", key="color_marker")
        color_enzyme = st.color_picker("Enzyme", key="color_enzyme")
        color_promoter = st.color_picker("Promoter", key="color_promoter")
        color_probe = st.color_picker("Probe", key="color_probe")
        color_line = st.color_picker("Backbone", key="color_line")

    with st.expander("📏 요소 두께 (mm)"):
        box_h_mm = st.slider("유전자 박스 높이", 1.0, 10.0, step=0.5, key="box_h_mm")
        probe_h_mm = st.slider("Probe 높이", 1.0, 10.0, step=0.5, key="probe_h_mm")
        enzyme_tick_mm = st.slider("Enzyme 눈금 길이", 1.0, 12.0, step=0.5, key="enzyme_tick_mm")
        primer_stem_mm = st.slider("Primer 세로선 길이", 1.0, 15.0, step=0.5, key="primer_stem_mm")
        arrow_mm = st.slider("Primer 화살표 길이", 1.0, 10.0, step=0.5, key="arrow_mm")
        exon_head_ratio = st.slider("Exon 화살촉 폭 (박스 대비)", 1.0, 3.0, step=0.1,
                                    key="exon_head_ratio")
        track_gap_mm = st.slider("WT–Mutant 최소 간격", 2.0, 30.0, step=1.0, key="track_gap_mm")
        line_width = st.slider("Backbone 두께 (pt)", 0.5, 5.0, step=0.5, key="line_width")

    with st.expander("✂️ 양 끝 자르기 (WT/Mutant 동시 적용)"):
        st.caption("자동 계산된 범위에서 양쪽 끝을 추가로 잘라냅니다. 0이면 비활성. 단위는 bp.")
        trim_left = st.number_input("⬅️ 왼쪽 (bp)", min_value=0, step=100, key="trim_left")
        trim_right = st.number_input("➡️ 오른쪽 (bp)", min_value=0, step=100, key="trim_right")

    with st.expander("📏 Scale bar"):
        show_scalebar = st.checkbox("Scale bar 표시", key="show_scalebar")
        scalebar_bp = st.number_input("길이 (bp)", min_value=100, step=100, key="scalebar_bp")
        scalebar_color = st.color_picker("색상", key="scalebar_color")
        scalebar_lw = st.slider("두께", 1.0, 6.0, step=0.5, key="scalebar_lw")
        st.caption("1000 bp 이상은 자동으로 kb로 표기됩니다.")

    with st.expander("🔧 겹침 처리 & 기타", expanded=True):
        auto_stagger = st.checkbox("라벨 겹침 자동 해소", key="auto_stagger")
        st.caption(
            "먼저 한 줄 안에서 좌우로만 밀어 겹침을 없애고(원래 위치에서 5 mm 이상 밀려야 하는 "
            "라벨만 윗단으로 올림), 많이 밀린 라벨에는 가는 지시선을 그립니다."
        )
        dedup_enzyme = st.checkbox("같은 효소는 첫 사이트에만 라벨", key="dedup_enzyme")
        st.caption("눈금은 모두 그리고 이름만 한 번 표기하는, 제한효소 지도의 표준 방식입니다.")
        exon_arrow = st.checkbox("마지막 exon을 화살표로", key="exon_arrow")
        st.caption("유전자 방향(strand)을 보고 끝쪽 exon의 끝을 뾰족하게 그립니다.")
        show_hr = st.checkbox("HR 점선 표시", key="show_hr")

    with st.expander("🧠 자동 인식 규칙"):
        st.caption(
            "이름에 아래 단어가 들어가면 해당 종류로 분류합니다. 쉼표로 구분하며 대소문자는 무시합니다. "
            "여기에 없어도 버려지지 않고, 아래 길이 기준으로 추정한 뒤 '요소 편집' 탭에서 바로 바꿀 수 있습니다."
        )
        rule_marker = st.text_input("Marker 단어", key="rule_marker")
        rule_promoter = st.text_input("Promoter 단어", key="rule_promoter")
        rule_probe = st.text_input("Probe 단어", key="rule_probe")
        rule_primer = st.text_input("Primer 단어", key="rule_primer")
        st.caption("이름으로 못 찾은 feature는 길이로 추정합니다.")
        rule_enzyme_max = st.number_input("이 길이 이하는 Enzyme (bp)", 1, 100, step=1, key="rule_enzyme_max")
        rule_primer_max = st.number_input("이 길이 이하는 Primer (bp)", 10, 500, step=5, key="rule_primer_max")
        hide_unnamed = st.checkbox("이름 없는 feature는 처음에 꺼두기", key="hide_unnamed")

    with st.expander("💾 프리셋 (내 기본값)"):
        current_preset = {k: st.session_state[k] for k in PRESET_DEFAULTS}
        st.download_button(
            "⬇️ 현재 설정을 JSON으로 저장",
            data=json.dumps(current_preset, indent=2, ensure_ascii=False),
            file_name="figure_preset.json",
            mime="application/json",
            use_container_width=True,
        )
        uploaded_preset = st.file_uploader("⬆️ 프리셋 불러오기", type=["json"], key="preset_uploader")
        if uploaded_preset is not None:
            sig = uploaded_preset.name + str(uploaded_preset.size)
            if st.session_state.get("_last_preset_sig") != sig:
                try:
                    loaded = json.loads(uploaded_preset.getvalue().decode("utf-8"))
                    for k, v in loaded.items():
                        if k in PRESET_DEFAULTS:
                            st.session_state[k] = v
                    st.session_state["_last_preset_sig"] = sig
                    st.success("프리셋을 적용했습니다.")
                    st.rerun()
                except Exception as e:
                    st.error(f"프리셋 파일을 읽지 못했습니다: {e}")
        if st.button("🔄 공장 초기값으로 리셋", use_container_width=True):
            for k, v in PRESET_DEFAULTS.items():
                st.session_state[k] = v
            st.session_state.pop("_last_preset_sig", None)
            st.rerun()

    st.markdown("---")
    st.caption(
        "⚠️ **As-is 제공** · 결과물 검증 책임은 사용자에게 있습니다.\n\n"
        "🔒 업로드 파일은 서버에 저장되지 않고 세션 종료 시 파기됩니다.\n\n"
        "📜 Noncommercial Academic Use / 비영리 학술 사용"
    )
    st.caption("👨‍💻 Developed by Yujin Lee | 💡 버그 제보 및 기능 건의 환영!")


# ============================================================
# 텍스트 서식
# ============================================================
def _format_text_legacy(text):
    """***굵은기울임*** / **굵게** / *기울임* / ^위첨자^ / ~아래첨자~ 를 mathtext로 변환."""
    if not text:
        return ""
    text = str(text)
    if "$" in text:
        return text
    if not re.search(r"\*\*\*|\*\*|\*|\^.+?\^|~.+?~", text):
        return text

    s = text
    prev = None
    while prev != s:
        prev = s
        s = re.sub(r"\*\*\*([^*]+?)\*\*\*", r"\\mathbf{\\mathit{\1}}", s)
        s = re.sub(r"\*\*([^*]+?)\*\*", r"\\mathbf{\1}", s)
        s = re.sub(r"(?<!\*)\*([^*]+?)\*(?!\*)", r"\\mathit{\1}", s)
        s = re.sub(r"\^([^\^]+?)\^", r"^{\1}", s)
        s = re.sub(r"~([^~]+?)~", r"_{\1}", s)
    s = s.replace(" ", r"\ ")
    return f"${s}$"


def _parse_text_markup(text, bold=False, italic=False, script=""):
    """수동 서식과 아래/위첨자를 문자별 속성으로 읽는다."""
    chars, i = [], 0
    while i < len(text):
        token = next((mark for mark in ("***", "**", "*", "~", "^") if text.startswith(mark, i)), None)
        if token:
            end = text.find(token, i + len(token))
            if end > i + len(token):
                b, it, sc = bold, italic, script
                if token == "***": b, it = True, True
                elif token == "**": b = True
                elif token == "*": it = True
                elif token == "~": sc = "sub"
                elif token == "^": sc = "sup"
                chars.extend(_parse_text_markup(text[i + len(token):end], b, it, sc))
                i = end + len(token)
                continue
        chars.append({"char": text[i], "bold": bold, "italic": italic, "script": script})
        i += 1
    return chars


def _mathtext_escape(text):
    escapes = {" ": r"\ ", "_": r"\_", "$": r"\$", "%": r"\%",
               "{": r"\{", "}": r"\}", "\\": r"\backslash ", "#": r"\#", "&": r"\&"}
    return "".join(escapes.get(char, char) for char in text)


def format_text(text, styles=None):
    """선택한 단어의 기울임을 덮어쓰며 기존 굵기·아래/위첨자는 보존한다."""
    if not styles or not text or "$" in str(text):
        return _format_text_legacy(text)
    chars = _parse_text_markup(str(text))
    plain = "".join(char["char"] for char in chars)
    for rule in styles:
        word = str(rule.get("단어", ""))
        if not word:
            continue
        mode = rule.get("일치 방식", "부분 일치")
        pattern = re.escape(word)
        if mode == "단어 일치": pattern = r"(?<!\w)" + pattern + r"(?!\w)"
        elif mode == "전체 일치": pattern = r"\A" + pattern + r"\Z"
        for match in re.finditer(pattern, plain, flags=re.IGNORECASE):
            for index in range(match.start(), match.end()):
                chars[index]["italic"] = bool(rule["italic"])
    if not any(char["bold"] or char["italic"] or char["script"] for char in chars):
        return plain
    groups, i = [], 0
    while i < len(chars):
        script = chars[i]["script"]
        end = i + 1
        while end < len(chars) and chars[end]["script"] == script: end += 1
        chunks, j = [], i
        while j < end:
            bold, italic = chars[j]["bold"], chars[j]["italic"]
            finish = j + 1
            while finish < end and (chars[finish]["bold"], chars[finish]["italic"]) == (bold, italic): finish += 1
            content = _mathtext_escape("".join(char["char"] for char in chars[j:finish]))
            if bold and italic: styled = r"\mathbf{\mathit{" + content + "}}"
            elif bold: styled = r"\mathbf{" + content + "}"
            elif italic: styled = r"\mathit{" + content + "}"
            else: styled = r"\mathrm{" + content + "}"
            chunks.append(styled)
            j = finish
        group = "".join(chunks)
        if script: group = ("_" if script == "sub" else "^") + "{" + group + "}"
        if not groups and script: groups.append("{}")
        groups.append(group)
        i = end
    return "$" + "".join(groups) + "$"


def rules_for_design(styles, design_id):
    return [rule for rule in (styles or []) if rule["디자인 ID"] in ("*", design_id)]


def styled_design_text(text, styles, target):
    applicable = [rule for rule in (styles or []) if rule["적용 대상"] in ("전체", target)]
    return format_text(text, applicable)


def plain_len(text):
    """서식 기호를 뺀 실제 글자 수 (폭 추정용)."""
    return max(len(re.sub(r"[\*\^~$]", "", str(text))), 1)


# ============================================================
# GenBank 파서
# ============================================================
def split_rule(text):
    return tuple(w.strip().upper() for w in str(text).split(",") if w.strip())


def classify(name, ftype, length, rules, all_enzymes):
    """이름 → 길이 순으로 추정. 어떤 경우에도 None을 돌려주지 않는다."""
    upper = name.upper()
    nospace = upper.replace(" ", "")
    if any(k in upper for k in rules["probe"]):
        return "Probe"
    if ftype == "restriction_site" or upper in all_enzymes or nospace in all_enzymes:
        return "Enzyme"
    if any(k in upper for k in rules["marker"]):
        return "Marker"
    if ftype == "promoter" or any(k in upper for k in rules["promoter"]):
        return "Promoter"
    if ftype in ("primer_bind", "misc_binding") or any(k in upper for k in rules["primer"]):
        return "Primer"
    if ftype in ("cds", "exon", "gene", "mrna", "trna"):
        return "Exon"
    if length <= rules["enzyme_max"]:
        return "Enzyme"
    if length <= rules["primer_max"]:
        return "Primer"
    return "Exon"


@st.cache_data
def parse_to_dataframe(file_content, is_mutant, rules_key):
    """GenBank의 모든 feature를 읽어 온다. source만 빼고 어떤 것도 버리지 않는다."""
    rules = {
        "marker": rules_key[0], "promoter": rules_key[1],
        "probe": rules_key[2], "primer": rules_key[3],
        "enzyme_max": rules_key[4], "primer_max": rules_key[5],
    }
    hide_unnamed = rules_key[6]

    record = SeqIO.read(StringIO(file_content), "genbank")
    try:
        all_enzymes = {str(e).upper() for e in Restriction.AllEnzymes}
    except Exception:
        all_enzymes = {"ECORV", "PSTI", "BAMHI", "ECORI", "HINDIII", "XHOI", "NOTI", "PPUMI"}

    rows, notes = [], []
    for f in record.features:
        ftype = f.type.lower()
        if ftype == "source":
            continue

        name, named = "", True
        for key in ["label", "note", "gene", "name", "ApEinfo_label", "product", "standard_name"]:
            if key in f.qualifiers and str(f.qualifiers[key][0]).strip():
                name = str(f.qualifiers[key][0]).strip()
                break
        if not name:
            name, named = f.type, False

        start, end = int(f.location.start), int(f.location.end)
        strand = f.location.strand if f.location.strand is not None else 1
        category = classify(name, ftype, end - start, rules, all_enzymes)

        if not named:
            notes.append(f"이름 없음 · {f.type} {start}–{end} → {category}")
        elif category in ("Enzyme", "Primer", "Exon") and ftype not in (
            "restriction_site", "primer_bind", "misc_binding", "cds", "exon", "gene", "mrna", "trna"
        ) and not any(k in name.upper() for k in rules["primer"] + rules["marker"]) \
                and name.upper() not in all_enzymes and name.upper().replace(" ", "") not in all_enzymes:
            notes.append(f"{name} ({end - start} bp) → 길이로 {category} 추정")

        if category == "Enzyme":
            up_default = False            # 효소는 선 아래 한 줄
        elif category == "Primer":
            up_default = True             # 프라이머는 선 위 한 줄
        elif category == "Probe":
            up_default = not is_mutant    # WT는 위, Mutant는 아래
        else:
            up_default = True

        show = named or not hide_unnamed
        if name.upper() in ("QLP", "QRP"):
            show = False

        rows.append({
            "표시": show,
            "이름": name,
            "종류": category,
            "화살표": "자동",
            "길이_배수": 1.0,
            "Y_띄우기_mm": 0.0,
            "위로향함": up_default,
            "시작": start, "종료": end, "방향": strand,
        })

    df = pd.DataFrame(rows, columns=["표시", "이름", "종류", "화살표", "길이_배수",
                                     "Y_띄우기_mm", "위로향함", "시작", "종료", "방향"])
    if not df.empty:
        df = df.sort_values("시작").reset_index(drop=True)
    return df, len(record), notes


# ============================================================
# 배치 계산 헬퍼
# ============================================================
def text_bp_width(text, fs, bp_per_inch):
    """라벨의 가로 폭을 bp 단위로 추정."""
    return plain_len(text) * fs * 0.62 / 72.0 * bp_per_inch


def text_height(fs):
    """라벨의 세로 높이를 inch 단위로 추정."""
    return fs * 1.35 / 72.0


def pack_row(items, gap_bp):
    """한 줄 안에서 라벨이 겹치지 않도록 가로로만 밀어낸다.
    원래 위치에서의 이동량 제곱합을 최소화하는 배치(등위회귀, PAVA)."""
    n = len(items)
    if n == 0:
        return {}
    W = [w + gap_bp for _, _, w in items]
    prefix, acc = [], 0.0
    for w in W:
        prefix.append(acc)
        acc += w
    u = [items[i][1] - W[i] / 2 - prefix[i] for i in range(n)]

    vals, cnts = [], []
    for x in u:
        vals.append(x)
        cnts.append(1)
        while len(vals) > 1 and vals[-2] > vals[-1]:
            v2, c2 = vals.pop(), cnts.pop()
            v1, c1 = vals.pop(), cnts.pop()
            vals.append((v1 * c1 + v2 * c2) / (c1 + c2))
            cnts.append(c1 + c2)
    v = []
    for val, c in zip(vals, cnts):
        v.extend([val] * c)

    return {items[i][0]: v[i] + prefix[i] + W[i] / 2 for i in range(n)}


def place_labels(items, fs, bp_per_inch, gap_bp, max_shift, enabled, max_tiers=3):
    """items: [(idx, x_ideal, label)] → {idx: (tier, x_label)}.
    먼저 한 줄에 다 넣어보고, 너무 멀리 밀려나는 것만 윗단으로 올린다."""
    if not items:
        return {}
    if not enabled:
        return {i: (0, x) for i, x, _ in items}

    result = {}
    remaining = sorted(items, key=lambda t: t[1])
    for tier in range(max_tiers):
        sized = [(i, x, text_bp_width(l, fs, bp_per_inch)) for i, x, l in remaining]
        packed = pack_row(sized, gap_bp)
        last = tier == max_tiers - 1
        pushed = []
        for i, x, l in remaining:
            if last or abs(packed[i] - x) <= max_shift:
                result[i] = (tier, packed[i])
            else:
                pushed.append((i, x, l))
        if not pushed or len(pushed) == len(remaining):
            for i, x, l in pushed:
                result[i] = (tier, packed[i])
            break
        remaining = pushed
    return result


def enzyme_label_flags(df, dedup):
    """같은 이름의 효소는 가장 왼쪽 사이트에만 라벨을 붙인다."""
    flags = {}
    seen = set()
    sub = df[(df["표시"] == True) & (df["종류"] == "Enzyme")].sort_values("시작")
    for idx, row in sub.iterrows():
        key = str(row["이름"]).upper()
        if dedup and key in seen:
            flags[idx] = False
        else:
            flags[idx] = True
            seen.add(key)
    return flags


def plan_track(df, cfg):
    """트랙 한 줄의 모든 요소를 배치. 좌표는 트랙 선으로부터의 '거리(inch)'로 저장."""
    plan = {}
    extent = {1: cfg["box_h"] / 2.0, -1: cfg["box_h"] / 2.0}
    vis = df[df["표시"] == True]
    flags = enzyme_label_flags(df, cfg["dedup_enzyme"])
    gap_bp = 1.0 * MM * cfg["bp_per_inch"]
    max_shift = 2.5 * MM * cfg["bp_per_inch"]   # 이보다 더 밀려야 하면 윗단으로
    arrow_dx = cfg["arrow_mm"] * MM * cfg["bp_per_inch"]

    for side in (1, -1):
        base = 0.0

        # --- Enzyme ---
        rows = [(i, r) for i, r in vis.iterrows()
                if r["종류"] == "Enzyme" and (1 if r["위로향함"] else -1) == side]
        labeled = [(i, float(r["시작"]), str(r["이름"])) for i, r in rows if flags.get(i, True)]
        placed = place_labels(labeled, cfg["enzyme_font_size"], cfg["bp_per_inch"],
                              gap_bp, max_shift, cfg["auto_stagger"])
        th = text_height(cfg["enzyme_font_size"])
        top = base
        for i, r in rows:
            x0 = float(r["시작"])
            tier, xlab = placed.get(i, (0, x0))
            y0 = base + float(r["Y_띄우기_mm"]) * MM
            y1 = y0 + cfg["enzyme_tick"] * float(r["길이_배수"]) + tier * (th + 0.6 * MM)
            has_label = flags.get(i, True)
            plan[i] = {"kind": "Enzyme", "side": side, "y0": y0, "y1": y1,
                       "label": has_label, "xlab": xlab, "xtip": x0}
            top = max(top, y1 + (th + 1.2 * MM if has_label else 0.0))
        if rows:
            base = top + 1.0 * MM

        # --- Primer ---
        rows = [(i, r) for i, r in vis.iterrows()
                if r["종류"] == "Primer" and (1 if r["위로향함"] else -1) == side]
        labeled = [(i, float(r["시작"]) + arrow_dx / 2 * (1 if r["방향"] == 1 else -1),
                    str(r["이름"])) for i, r in rows]
        placed = place_labels(labeled, cfg["primer_font_size"], cfg["bp_per_inch"],
                              gap_bp, max_shift, cfg["auto_stagger"])
        th = text_height(cfg["primer_font_size"])
        top = base
        for i, r in rows:
            ideal = float(r["시작"]) + arrow_dx / 2 * (1 if r["방향"] == 1 else -1)
            tier, xlab = placed.get(i, (0, ideal))
            y0 = base + float(r["Y_띄우기_mm"]) * MM
            y1 = y0 + cfg["primer_stem"] * float(r["길이_배수"]) + tier * (th + 0.6 * MM)
            plan[i] = {"kind": "Primer", "side": side, "y0": y0, "y1": y1,
                       "label": True, "xlab": xlab, "xtip": ideal, "dx": arrow_dx}
            top = max(top, y1 + th + 1.2 * MM)
        if rows:
            base = top + 1.0 * MM

        # --- Probe (해당 면의 가장 바깥, 같은 면끼리 높이·두께·글씨 크기 완전 동일) ---
        rows = [(i, r) for i, r in vis.iterrows()
                if r["종류"] == "Probe" and (1 if r["위로향함"] else -1) == side]
        if rows:
            fs = cfg["marker_font_size"]
            th = text_height(fs)
            # 모든 Probe 이름은 박스 바로 아래: 위쪽 Probe는 먼저 라벨 공간 확보.
            y0 = base + 1.2 * MM + (th + 0.8 * MM if side == 1 else 0.0)
            y1 = y0 + cfg["probe_h"]
            for i, r in rows:
                plan[i] = {"kind": "Probe", "side": side, "y0": y0, "y1": y1,
                           "label": True, "fontsize": fs, "label_outside": True}
            base = y1 + 1.0 * MM + (th + 0.8 * MM if side == -1 else 0.0)

        extent[side] = max(extent[side], base)

    # --- 박스류: 라벨이 박스보다 길면 6pt까지 축소, 그래도 넘치면 스택 바깥으로 ---
    exon_rows = [(i, r) for i, r in vis.iterrows() if r["종류"] == "Exon"]
    arrow_idx, arrow_dir = None, 1
    if cfg.get("exon_arrow") and exon_rows:
        strands = [int(r["방향"]) for _, r in exon_rows]
        arrow_dir = 1 if sum(1 for x in strands if x >= 0) >= len(strands) / 2 else -1
        key = (lambda t: float(t[1]["종료"])) if arrow_dir == 1 else (lambda t: -float(t[1]["시작"]))
        arrow_idx = max(exon_rows, key=key)[0]
    any_arrow = False

    overflow, long_names = [], []
    for i, r in vis.iterrows():
        if r["종류"] not in ("Exon", "Marker", "Promoter"):
            continue
        w = float(r["종료"]) - float(r["시작"])

        choice = str(r.get("화살표", "자동"))
        if choice == "→":
            arrow, adir = True, 1
        elif choice == "←":
            arrow, adir = True, -1
        elif choice == "없음":
            arrow, adir = False, 1
        else:                                   # 자동
            if r["종류"] == "Marker":
                # NAT/NEO/HYG 등 마커는 각 feature의 strand 방향으로 화살표 표시
                arrow, adir = True, (-1 if r["방향"] == -1 else 1)
            else:
                arrow, adir = (i == arrow_idx), arrow_dir
        any_arrow = any_arrow or arrow

        # 화살촉 부분을 뺀 자루 폭 기준으로 라벨을 앉힌다
        tip = min(cfg["box_h"] * cfg.get("head_ratio", 1.8) * cfg["bp_per_inch"],
                  w * 0.45) if arrow else 0.0
        shaft = w - tip
        fs = cfg["marker_font_size"]
        outside = False
        if r["종류"] in ("Marker", "Promoter"):
            need = text_bp_width(r["이름"], fs, cfg["bp_per_inch"])
            if need > shaft * 0.92:
                fs = max(6.0, fs * (shaft * 0.92) / max(need, 1e-9))
                if text_bp_width(r["이름"], fs, cfg["bp_per_inch"]) > shaft * 0.92:
                    outside = True
                    fs = cfg["marker_font_size"]
                    overflow.append((i, float(r["시작"]) + w / 2, str(r["이름"])))
                    long_names.append(str(r["이름"]))
        plan[i] = {"kind": r["종류"], "side": cfg["out_side"], "fontsize": fs,
                   "outside": outside, "arrow": arrow, "arrow_dir": adir,
                   "tip": tip}

    if any_arrow:   # 화살촉이 박스보다 두꺼우므로 위아래 여유 확보
        head_half = cfg["box_h"] * cfg.get("head_ratio", 1.8) / 2.0
        extent[1] = max(extent[1], head_half)
        extent[-1] = max(extent[-1], head_half)

    if overflow:
        th = text_height(cfg["marker_font_size"])
        placed = place_labels(overflow, cfg["marker_font_size"], cfg["bp_per_inch"],
                              gap_bp, max_shift, True, max_tiers=2)
        max_tier = max((tier for tier, _ in placed.values()), default=0)
        if cfg.get("outside_below", False):
            # Mutant: 긴 이름은 해당 박스 바로 아래에 표시한다.
            side_out = -1
            box_bottom = max(
                cfg["box_h"] * (cfg.get("head_ratio", 1.8) if plan[i]["arrow"] else 1.0) / 2
                for i, _, _ in overflow
            )
            y_out = box_bottom + 0.6 * MM
            label_bottom = y_out + max_tier * (th + 0.6 * MM) + th
            # 효소 눈금은 시퀀스 선에 붙여 둔다. Probe만 필요할 때 간격을 확보한다.
            bottom_probes = [geo for geo in plan.values()
                             if geo["side"] == -1 and geo["kind"] == "Probe"]
            probe_shift = max(0.0, label_bottom + 1.2 * MM
                              - min(geo["y0"] for geo in bottom_probes)) if bottom_probes else 0.0
            for geo in bottom_probes:
                geo["y0"] += probe_shift
                geo["y1"] += probe_shift
            extent[-1] = max(extent[-1] + probe_shift, label_bottom)
        else:
            side_out = 1 if extent[1] <= extent[-1] else -1
            y_out = extent[side_out] + 1.2 * MM
            extent[side_out] = y_out + max_tier * (th + 0.6 * MM) + th
        for i, x, name in overflow:
            tier, xlab = placed.get(i, (0, x))
            plan[i].update(side=side_out, y_out=y_out + tier * (th + 0.6 * MM), xlab=xlab)

    return plan, extent, long_names


def render_gene_map(wt_df, mut_df, wt_len, mut_len, hr_edited, settings,
                    wt_label="WT allele", mut_label="Mutant allele", text_styles=None):
    """단일/일괄 작업에서 같은 레이아웃으로 Figure를 생성한다."""
    trim_left = settings["trim_left"]
    trim_right = settings["trim_right"]
    box_h_mm = settings["box_h_mm"]
    probe_h_mm = settings["probe_h_mm"]
    enzyme_tick_mm = settings["enzyme_tick_mm"]
    primer_stem_mm = settings["primer_stem_mm"]
    label_area_mm = settings["label_area_mm"]
    fig_width = settings["fig_width"]
    marker_font_size = settings["marker_font_size"]
    enzyme_font_size = settings["enzyme_font_size"]
    primer_font_size = settings["primer_font_size"]
    arrow_mm = settings["arrow_mm"]
    dedup_enzyme = settings["dedup_enzyme"]
    auto_stagger = settings["auto_stagger"]
    exon_arrow = settings["exon_arrow"]
    exon_head_ratio = settings["exon_head_ratio"]
    track_gap_mm = settings["track_gap_mm"]
    show_scalebar = settings["show_scalebar"]
    scalebar_lw = settings["scalebar_lw"]
    label_font_size = settings["label_font_size"]
    color_line = settings["color_line"]
    line_width = settings["line_width"]
    text_color = settings["text_color"]
    color_enzyme = settings["color_enzyme"]
    color_exon = settings["color_exon"]
    color_marker = settings["color_marker"]
    color_promoter = settings["color_promoter"]
    color_probe = settings["color_probe"]
    show_hr = settings["show_hr"]
    scalebar_bp = settings["scalebar_bp"]
    scalebar_color = settings["scalebar_color"]

    # ---------- 표시 범위 ----------
    all_visible = pd.concat([wt_df[wt_df["표시"] == True], mut_df[mut_df["표시"] == True]])
    if not all_visible.empty:
        data_start = int(all_visible["시작"].min())
        data_end = int(all_visible["종료"].max())
    else:
        data_start, data_end = 0, max(wt_len, mut_len)

    pad = max(int((data_end - data_start) * 0.02), 50)
    final_start = max(0, data_start - pad) + trim_left
    final_end = (data_end + pad) - trim_right
    if final_end <= final_start:
        final_end = final_start + 100
    span = final_end - final_start

    # ---------- 좌표계: x는 bp, y는 inch (축이 figure 전체를 채움) ----------
    box_h = box_h_mm * MM
    probe_h = probe_h_mm * MM
    enzyme_tick = enzyme_tick_mm * MM
    primer_stem = primer_stem_mm * MM
    label_area = label_area_mm * MM
    right_pad = 0.12

    draw_width = max(fig_width - label_area - right_pad, 1.0)
    bp_per_inch = span / draw_width

    cfg_common = dict(
        span=span, bp_per_inch=bp_per_inch,
        box_h=box_h, probe_h=probe_h,
        enzyme_tick=enzyme_tick, primer_stem=primer_stem,
        enzyme_font_size=enzyme_font_size, primer_font_size=primer_font_size,
        marker_font_size=marker_font_size, arrow_mm=arrow_mm,
        dedup_enzyme=dedup_enzyme, auto_stagger=auto_stagger, exon_arrow=exon_arrow,
        head_ratio=exon_head_ratio,
    )

    wt_plan, wt_ext, wt_long = plan_track(wt_df, dict(cfg_common, out_side=1))
    mut_plan, mut_ext, mut_long = plan_track(mut_df, dict(cfg_common, out_side=-1, outside_below=True))

    gap = wt_ext[-1] + mut_ext[1] + track_gap_mm * MM
    y_wt = 0.0
    y_mut = -gap

    y_top = y_wt + wt_ext[1]
    y_bottom = y_mut - mut_ext[-1]

    if show_scalebar:
        y_top += 4 * MM + scalebar_lw / 72.0 + text_height(label_font_size)

    # WT/Mutant 라벨이 트랙보다 크면 위아래 여유 확보
    half_label = text_height(label_font_size) / 2.0
    y_top = max(y_top, y_wt + half_label)
    y_bottom = min(y_bottom, y_mut - half_label)

    margin = 2 * MM
    y_lo, y_hi = y_bottom - margin, y_top + margin
    fig_height = max(y_hi - y_lo, 0.8)

    fig = plt.figure(figsize=(fig_width, fig_height))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(final_start - label_area * bp_per_inch, final_end + right_pad * bp_per_inch)
    ax.set_ylim(y_lo, y_hi)
    ax.axis("off")

    Z_LINE, Z_HR, Z_BOX, Z_MARK = 1, 2, 5, 8

    for y in (y_wt, y_mut):
        ax.plot([final_start, final_end], [y, y], color=color_line, lw=line_width,
                solid_capstyle="butt", zorder=Z_LINE)

    ax.text(final_start - 2 * MM * bp_per_inch, y_wt, styled_design_text(wt_label, text_styles, "WT 라벨"),
            ha="right", va="center", color=text_color, fontsize=label_font_size)
    ax.text(final_start - 2 * MM * bp_per_inch, y_mut, styled_design_text(mut_label, text_styles, "Mutant 라벨"),
            ha="right", va="center", color=text_color, fontsize=label_font_size)


    def draw_track(track_y, df, plan, target):
        dropped = []
        for idx, row in df[df["표시"] == True].iterrows():
            geo = plan.get(idx)
            if geo is None:
                continue
            s, e = float(row["시작"]), float(row["종료"])
            if e <= final_start or s >= final_end:
                dropped.append(str(row["이름"]))
                continue
            kind = geo["kind"]
            name = styled_design_text(row["이름"], text_styles, target)

            if kind in ("Enzyme", "Primer"):
                if s < final_start or s > final_end:
                    dropped.append(str(row["이름"]))
                    continue
                side = geo["side"]
                y0 = track_y + side * geo["y0"]
                y1 = track_y + side * geo["y1"]
                va = "bottom" if side == 1 else "top"

                xlab = geo.get("xlab", s)
                xtip = geo.get("xtip", s)
                y_text = y1 + side * 1.0 * MM

                if kind == "Enzyme":
                    ax.plot([s, s], [y0, y1], color=color_enzyme, lw=1.0,
                            solid_capstyle="butt", zorder=Z_BOX - 1)
                    if geo["label"]:
                        ax.text(xlab, y_text, name, color=color_enzyme, ha="center", va=va,
                                fontsize=enzyme_font_size, zorder=Z_MARK)
                else:
                    dx = geo.get("dx", 2.5 * MM * bp_per_inch) * (1 if row["방향"] == 1 else -1)
                    ax.plot([s, s], [y0, y1], color=color_line, lw=0.9,
                            solid_capstyle="butt", zorder=Z_BOX - 1)
                    ax.annotate("", xy=(s + dx, y1), xytext=(s, y1),
                                arrowprops=dict(arrowstyle="-|>", color=color_line,
                                                lw=0.9, mutation_scale=7,
                                                shrinkA=0, shrinkB=0),
                                zorder=Z_BOX - 1)
                    ax.text(xlab, y_text, name, ha="center", va=va,
                            color=text_color, fontsize=primer_font_size, zorder=Z_MARK)
                continue

            # 박스류 — 표시 범위 경계에 맞춰 자르되, 최소 폭은 보장해서 사라지지 않게
            sx, ex = max(s, final_start), min(e, final_end)
            w = max(ex - sx, 0.35 * MM * bp_per_inch)

            if kind == "Probe":
                side = geo["side"]
                y0 = track_y + side * geo["y0"]
                y1 = track_y + side * geo["y1"]
                lo, hi = min(y0, y1), max(y0, y1)
                fs = geo.get("fontsize", marker_font_size)
                ax.add_patch(patches.Rectangle((sx, lo), w, hi - lo,
                                               facecolor=color_probe, edgecolor="none",
                                               zorder=Z_BOX))
                # 박스 폭이나 트랙의 위/아래 위치와 관계없이 이름은 항상 박스 아래.
                ax.text(sx + w / 2, lo - 0.8 * MM, name, ha="center", va="top",
                        fontsize=fs, color=text_color, zorder=Z_MARK)

            elif kind in ("Exon", "Marker", "Promoter"):
                face = {"Exon": color_exon, "Marker": color_marker, "Promoter": color_promoter}[kind]
                d = geo.get("arrow_dir", 1)
                # 화살표 끝이 표시 범위 안에 있으면 화살촉을 그린다
                tip_ok = geo.get("arrow") and (e <= final_end if d == 1 else s >= final_start)
                if tip_ok:
                    head_h = box_h * exon_head_ratio
                    tip = geo.get("tip") or min(head_h * bp_per_inch, w * 0.45)
                    lo, hi = track_y - box_h / 2, track_y + box_h / 2
                    hlo, hhi = track_y - head_h / 2, track_y + head_h / 2
                    if d == 1:
                        nx = sx + w - tip
                        pts = [(sx, lo), (nx, lo), (nx, hlo), (sx + w, track_y),
                               (nx, hhi), (nx, hi), (sx, hi)]
                    else:
                        nx = sx + tip
                        pts = [(sx + w, lo), (nx, lo), (nx, hlo), (sx, track_y),
                               (nx, hhi), (nx, hi), (sx + w, hi)]
                    ax.add_patch(patches.Polygon(pts, closed=True, facecolor=face,
                                                 edgecolor="none", zorder=Z_BOX))
                else:
                    ax.add_patch(patches.Rectangle((sx, track_y - box_h / 2), w, box_h,
                                                   facecolor=face, edgecolor="none", zorder=Z_BOX))
                if kind == "Exon":
                    continue
                fs = geo.get("fontsize", marker_font_size)
                if geo.get("outside"):
                    side = geo["side"]
                    ax.text(geo.get("xlab", sx + w / 2),
                            track_y + side * geo.get("y_out", box_h / 2 + 1 * MM), name,
                            ha="center", va="bottom" if side == 1 else "top",
                            fontsize=fs, color=text_color, zorder=Z_MARK)
                else:
                    inner = "#000000" if kind in ("Marker", "Promoter") else "#FFFFFF"
                    cx = sx + w / 2 - (d * geo.get("tip", 0.0) / 2 if tip_ok else 0.0)
                    ax.text(cx, track_y, name, ha="center", va="center",
                            fontsize=fs, color=inner, zorder=Z_MARK)
        return dropped


    dropped_wt = draw_track(y_wt, wt_df, wt_plan, "WT 내부")
    dropped_mut = draw_track(y_mut, mut_df, mut_plan, "Mutant 내부")

    # ---------- HR 점선 ----------
    if show_hr:
        for _, row in hr_edited.iterrows():
            wt_k = str(row.get("WT_이름", "")).strip()
            mut_k = str(row.get("Mutant_이름", "")).strip()
            if not wt_k or not mut_k or wt_k.lower() == "nan" or mut_k.lower() == "nan":
                continue
            wm = wt_df[(wt_df["표시"] == True) & (wt_df["이름"] == wt_k)]
            mm_ = mut_df[(mut_df["표시"] == True) & (mut_df["이름"] == mut_k)]
            if wm.empty or mm_.empty:
                continue
            wx, mx = float(wm.iloc[0]["시작"]), float(mm_.iloc[0]["시작"])
            if final_start <= wx <= final_end and final_start <= mx <= final_end:
                ax.plot([wx, mx], [y_wt, y_mut], "--", color=color_line,
                        lw=0.8, alpha=0.45, zorder=Z_HR)

    # ---------- Scale bar ----------
    if show_scalebar:
        sb_y = y_hi - margin - text_height(label_font_size)
        sb_x0, sb_x1 = final_start, final_start + scalebar_bp
        ax.plot([sb_x0, sb_x1], [sb_y, sb_y], color=scalebar_color, lw=scalebar_lw,
                solid_capstyle="butt", zorder=20)
        sb_label = f"{scalebar_bp / 1000:g} kb" if scalebar_bp >= 1000 else f"{scalebar_bp} bp"
        ax.text((sb_x0 + sb_x1) / 2, sb_y + 0.8 * MM, sb_label, ha="center", va="bottom",
                color=scalebar_color, fontsize=label_font_size, zorder=20)

    return fig, {"fig_height": fig_height, "dropped_wt": dropped_wt,
                 "dropped_mut": dropped_mut, "wt_long": wt_long, "mut_long": mut_long}


# ============================================================
# Excel 설정 양식 — 표준 XLSX 파일을 읽어 표시 설정만 적용
# ============================================================
EXCEL_DESIGN_COLUMNS = ["디자인 ID", "출력", "WT 파일", "Mutant 파일", "출력 이름", "WT 라벨", "Mutant 라벨"]
EXCEL_STYLE_COLUMNS = ["디자인 ID", "적용 대상", "단어", "기울임", "일치 방식"]
TEXT_STYLE_TARGETS = ["전체", "WT 라벨", "Mutant 라벨", "WT 내부", "Mutant 내부"]
TEXT_STYLE_MODES = ["부분 일치", "단어 일치", "전체 일치"]


def _cell_column(reference):
    column = 0
    for letter in re.match(r"[A-Z]+", reference.upper()).group(0):
        column = column * 26 + ord(letter) - ord("A") + 1
    return column - 1


def _read_xlsx_tables(raw):
    """shared strings / inline strings를 지원하는 읽기 전용 XLSX 판독기."""
    ns = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    rel_ns = "http://schemas.openxmlformats.org/package/2006/relationships"
    attr_ns = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    with zipfile.ZipFile(BytesIO(raw)) as archive:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        links = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        relationships = {item.attrib["Id"]: item.attrib["Target"] for item in links.findall("{" + rel_ns + "}Relationship")}
        strings = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            strings = ["".join(node.itertext()) for node in root.findall("s:si", ns)]
        tables = {}
        for sheet in workbook.findall("s:sheets/s:sheet", ns):
            name = sheet.attrib["name"]
            if name not in {"Designs", "TextStyles"}: continue
            target = relationships[sheet.attrib["{" + attr_ns + "}id"]]
            path = posixpath.normpath(target.lstrip("/") if target.startswith("/") else "xl/" + target)
            root = ET.fromstring(archive.read(path))
            rows = []
            for row in root.findall("s:sheetData/s:row", ns):
                values = {}
                for cell in row.findall("s:c", ns):
                    if cell.find("s:f", ns) is not None:
                        raise ValueError(f"{name} 시트에는 수식 대신 실제 텍스트를 입력해 주세요.")
                    kind = cell.attrib.get("t", "")
                    value = cell.find("s:v", ns)
                    value = value.text if value is not None else ""
                    if kind == "s": value = strings[int(value)]
                    elif kind == "inlineStr": value = "".join(node.text or "" for node in cell.findall("s:is//s:t", ns))
                    elif kind == "b": value = value == "1"
                    values[_cell_column(cell.attrib["r"])] = value
                rows.append(values)
            if not rows:
                tables[name] = pd.DataFrame()
                continue
            width = max(rows[0].keys(), default=-1) + 1
            headers = [str(rows[0].get(i, "")).strip() for i in range(width)]
            if len(headers) != len(set(headers)) or any(not header for header in headers):
                raise ValueError(f"{name} 1행의 열 이름을 양식 그대로 유지해 주세요.")
            tables[name] = pd.DataFrame([[row.get(i, "") for i in range(width)] for row in rows[1:]], columns=headers)
        return tables


def excel_text(value, default=""):
    if value is None or (not isinstance(value, (list, dict)) and pd.isna(value)): return default
    value = str(value).strip()
    return value or default


def excel_yes(value, default=True):
    text = excel_text(value).casefold()
    if not text: return default
    if text in {"예", "true", "yes", "y", "1"}: return True
    if text in {"아니오", "false", "no", "n", "0"}: return False
    raise ValueError(f"예 또는 아니오로 입력해 주세요: {value}")


@st.cache_data
def read_excel_configuration(raw):
    tables = _read_xlsx_tables(raw)
    if "Designs" not in tables: raise ValueError("Designs 시트가 필요합니다.")
    designs = tables["Designs"]
    missing = [name for name in EXCEL_DESIGN_COLUMNS if name not in designs.columns]
    if missing: raise ValueError("Designs에 필요한 열: " + ", ".join(missing))
    jobs, ids = [], set()
    for number, (_, row) in enumerate(designs.iterrows(), 2):
        if not any(excel_text(value) for value in row.values): continue
        design_id = excel_text(row["디자인 ID"])
        if not design_id or design_id == "*" or design_id in ids:
            raise ValueError(f"Designs {number}행의 디자인 ID는 고유한 값이어야 합니다.")
        ids.add(design_id)
        wt_file, mut_file = excel_text(row["WT 파일"]), excel_text(row["Mutant 파일"])
        if not wt_file or not mut_file: raise ValueError(f"Designs {number}행에 WT·Mutant 파일 이름을 입력해 주세요.")
        jobs.append({"디자인 ID": design_id, "출력": excel_yes(row["출력"]),
                     "WT 파일": wt_file, "Mutant 파일": mut_file,
                     "출력 이름": excel_text(row["출력 이름"], Path(mut_file).stem),
                     "WT 라벨": excel_text(row["WT 라벨"], "WT allele"),
                     "Mutant 라벨": excel_text(row["Mutant 라벨"], "Mutant allele")})
    if not jobs: raise ValueError("Designs에 디자인을 한 개 이상 입력해 주세요.")
    style_df = tables.get("TextStyles", pd.DataFrame(columns=EXCEL_STYLE_COLUMNS))
    missing = [name for name in EXCEL_STYLE_COLUMNS if name not in style_df.columns]
    if missing: raise ValueError("TextStyles에 필요한 열: " + ", ".join(missing))
    styles = []
    for number, (_, row) in enumerate(style_df.iterrows(), 2):
        word = excel_text(row["단어"])
        if not word: continue
        design_id = excel_text(row["디자인 ID"], "*")
        if design_id != "*" and design_id not in ids:
            raise ValueError(f"TextStyles {number}행의 디자인 ID가 Designs에 없습니다: {design_id}")
        target = excel_text(row["적용 대상"], "전체")
        mode = excel_text(row["일치 방식"], "부분 일치")
        if target not in TEXT_STYLE_TARGETS: raise ValueError(f"TextStyles {number}행의 적용 대상을 확인해 주세요.")
        if mode not in TEXT_STYLE_MODES: raise ValueError(f"TextStyles {number}행의 일치 방식을 확인해 주세요.")
        styles.append({"디자인 ID": design_id, "적용 대상": target, "단어": word,
                       "italic": excel_yes(row["기울임"]), "일치 방식": mode})
    return pd.DataFrame(jobs), styles



EXCEL_TEMPLATE_BASE64 = "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"

# ============================================================
# 일괄 작업 — 파일 짝짓기와 그림 ZIP 생성
# ============================================================
def decode_genbank_bytes(raw):
    for encoding in ("utf-8-sig", "latin-1"):
        try:
            return raw.decode(encoding)
        except UnicodeDecodeError:
            continue
    raise ValueError("GenBank 파일을 읽을 수 없습니다.")


def default_hr_pairs():
    return pd.DataFrame({"WT_이름": ["L1", "L2", "R1", "R2"],
                         "Mutant_이름": ["L2", "L1", "R2", "R1"]})


def file_pair_key(filename):
    """파일명만 사용해 후보를 추천한다. 모호한 짝은 선택하지 않는다."""
    stem = Path(filename).stem
    locus = re.search(r"CNAG[_ -]?(\d+)", stem, re.IGNORECASE)
    if locus:
        return "CNAG_" + locus.group(1)
    words = re.split(r"[\s_\-().]+", stem.upper())
    ignored = {"WT", "MUT", "MUTANT", "ALLELE", "WILD", "TYPE"}
    return "_".join(word for word in words if word and word not in ignored)


def suggest_batch_pairs(wt_files, mutant_files, wt_label="WT allele", mut_label="Mutant allele"):
    by_key = {}
    for uploaded in wt_files:
        by_key.setdefault(file_pair_key(uploaded.name), []).append(uploaded.name)
    rows = []
    for order, uploaded in enumerate(mutant_files, 1):
        candidates = by_key.get(file_pair_key(uploaded.name), [])
        candidate = candidates[0] if len(candidates) == 1 else ""
        rows.append({"디자인 ID": f"D{order:03d}", "출력": True, "WT 파일": candidate, "Mutant 파일": uploaded.name,
                     "출력 이름": Path(uploaded.name).stem,
                     "WT 라벨": wt_label, "Mutant 라벨": mut_label})
    return pd.DataFrame(rows, columns=EXCEL_DESIGN_COLUMNS)


def upload_lookup(uploaded_files):
    lookup = {}
    for uploaded in uploaded_files:
        if uploaded.name in lookup:
            raise ValueError(f"파일 이름이 중복됩니다: {uploaded.name}. 서로 다른 이름으로 바꿔 올려주세요.")
        lookup[uploaded.name] = uploaded
    return lookup


def safe_output_stem(value, fallback):
    name = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", str(value)).strip(" .")
    if not name or name.lower() == "nan":
        return fallback
    return name[:150]


def batch_signature(wt_files, mutant_files, pair_df, settings, formats, hr_pairs, rules_key, text_styles=None):
    digest = hashlib.sha256()
    for role, files in (("WT", wt_files), ("Mutant", mutant_files)):
        for uploaded in files:
            digest.update(role.encode())
            digest.update(uploaded.name.encode("utf-8"))
            digest.update(hashlib.sha256(uploaded.getvalue()).digest())
    digest.update(pair_df.to_json(force_ascii=False).encode("utf-8"))
    digest.update(hr_pairs.to_json(force_ascii=False).encode("utf-8"))
    digest.update(json.dumps([settings, formats, rules_key, text_styles], ensure_ascii=False, sort_keys=True).encode("utf-8"))
    return digest.hexdigest()


def build_batch_zip(pair_df, wt_files, mutant_files, rules_key, settings, formats,
                    hr_pairs=None, progress_callback=None, text_styles=None):
    wt_lookup, mut_lookup = upload_lookup(wt_files), upload_lookup(mutant_files)
    formats = [str(fmt).lower() for fmt in formats]
    if not formats or any(fmt not in {"png", "svg", "pdf"} for fmt in formats):
        raise ValueError("PNG, SVG, PDF 중 하나 이상의 형식을 선택해 주세요.")
    if hr_pairs is None:
        hr_pairs = default_hr_pairs()
    selected = pair_df[pair_df["출력"] == True]
    buffer, reports, used_names, parsed = BytesIO(), [], set(), {}
    success = 0

    def parse_file(uploaded, mutant):
        key = (mutant, uploaded.name)
        if key not in parsed:
            parsed[key] = parse_to_dataframe(decode_genbank_bytes(uploaded.getvalue()), mutant, rules_key)
        return parsed[key]

    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for order, (_, row) in enumerate(selected.iterrows(), 1):
            fig = None
            report = {"WT 파일": row.get("WT 파일", ""), "Mutant 파일": row.get("Mutant 파일", ""),
                      "출력 이름": row.get("출력 이름", ""), "상태": "실패", "내용": ""}
            try:
                wt_name, mut_name = str(row.get("WT 파일", "")), str(row.get("Mutant 파일", ""))
                if wt_name not in wt_lookup or mut_name not in mut_lookup:
                    raise ValueError("WT·Mutant 파일의 짝을 선택해 주세요.")
                wt_df, wt_len, wt_notes = parse_file(wt_lookup[wt_name], False)
                mut_df, mut_len, mut_notes = parse_file(mut_lookup[mut_name], True)
                fig, info = render_gene_map(
                    wt_df, mut_df, wt_len, mut_len, hr_pairs, settings,
                    str(row.get("WT 라벨", "WT allele")), str(row.get("Mutant 라벨", "Mutant allele")),
                    text_styles=rules_for_design(text_styles, str(row.get("디자인 ID", ""))))
                stem = safe_output_stem(row.get("출력 이름", ""), f"design_{order:03d}")
                root, suffix = stem, 2
                while stem.casefold() in used_names:
                    stem = f"{root}_{suffix}"
                    suffix += 1
                # 한 디자인의 모든 형식이 생성된 뒤 ZIP에 추가한다.
                entries = []
                for fmt in formats:
                    output = BytesIO()
                    kwargs = {"dpi": settings["dpi_setting"]} if fmt == "png" else {}
                    fig.savefig(output, format=fmt, bbox_inches="tight", pad_inches=0.05,
                                facecolor="white", **kwargs)
                    entries.append((f"{stem}/{stem}.{fmt}", output.getvalue()))
                for name, content in entries:
                    archive.writestr(name, content)
                used_names.add(stem.casefold())
                success += 1
                warnings = list(wt_notes) + list(mut_notes)
                clipped = info["dropped_wt"] + info["dropped_mut"]
                if clipped:
                    warnings.append("표시 범위 밖: " + ", ".join(clipped))
                if info["wt_long"] or info["mut_long"]:
                    warnings.append("박스 밖 라벨: " + ", ".join(info["wt_long"] + info["mut_long"]))
                report.update({"출력 이름": stem, "상태": "완료", "내용": " / ".join(warnings)})
            except Exception as error:
                report["내용"] = str(error)
            finally:
                if fig is not None:
                    plt.close(fig)
            reports.append(report)
            if progress_callback:
                progress_callback(order, len(selected))
        report_df = pd.DataFrame(reports, columns=["WT 파일", "Mutant 파일", "출력 이름", "상태", "내용"])
        archive.writestr("processing_report.csv", report_df.to_csv(index=False).encode("utf-8-sig"))
    return buffer.getvalue(), report_df, success


# ============================================================
# 메인
# ============================================================
tab_upload, tab_edit, tab_hr, tab_view, tab_free, tab_batch, tab_excel = st.tabs(
    ["📂 파일 업로드", "🛠️ 요소 편집", "🔗 HR 연결", "🖼️ 미리보기 & 다운로드", "✍️ 직접 편집", "📦 일괄 작업", "📋 엑셀 설정"]
)

with tab_upload:
    c1, c2 = st.columns(2)
    with c1:
        wt_file = st.file_uploader("WT GenBank (.gb)", type=["gb", "gbk"])
    with c2:
        mutant_file = st.file_uploader("Mutant GenBank (.gb)", type=["gb", "gbk"])

rules_key = (
    split_rule(rule_marker), split_rule(rule_promoter),
    split_rule(rule_probe), split_rule(rule_primer),
    int(rule_enzyme_max), int(rule_primer_max), bool(hide_unnamed),
)

excel_config, excel_styles = None, []
excel_upload = None
with tab_excel:
    st.caption("Designs에 파일별 이름을, TextStyles에 특정 단어의 기울임 여부를 입력하세요. 이 설정은 단일 그림과 일괄 작업 모두에 적용됩니다.")
    st.download_button("📥 엑셀 입력 양식 다운로드", data=base64.b64decode(EXCEL_TEMPLATE_BASE64),
        file_name="Gene_Map_Input_Template.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", key="excel_template")
    excel_upload = st.file_uploader("작성한 엑셀 설정 업로드", type=["xlsx"], key="excel_settings")
    if excel_upload:
        try:
            excel_config, excel_styles = read_excel_configuration(excel_upload.getvalue())
            st.success(f"{len(excel_config)}개 디자인과 {len(excel_styles)}개 서식 규칙을 불러왔습니다.")
            st.dataframe(excel_config, hide_index=True, use_container_width=True)
            if excel_styles:
                style_preview = pd.DataFrame(excel_styles).rename(columns={"italic": "기울임"})
                style_preview["기울임"] = style_preview["기울임"].map({True: "예", False: "아니오"})
                st.dataframe(style_preview, hide_index=True, use_container_width=True)
        except Exception as error:
            st.error(f"엑셀 설정을 읽지 못했습니다: {error}")
    st.caption("디자인 ID의 *는 모든 그림에 적용합니다. 아래쪽 규칙이 우선하며 기울임만 바뀌므로 굵기와 아래첨자는 유지됩니다. $...$ 직접 수식은 단어별 규칙에서 제외됩니다.")

# 일괄 작업은 단일 파일 업로드 없이도 사용할 수 있다.
with tab_batch:
    st.caption("WT·Mutant 파일을 여러 개 올리고 아래 표에서 짝과 라벨을 확인하세요. 현재 사이드바의 색상·크기·간격이 모든 그림에 적용됩니다.")
    batch_c1, batch_c2 = st.columns(2)
    with batch_c1:
        batch_wt = st.file_uploader("WT 파일 여러 개", type=["gb", "gbk"],
                                    accept_multiple_files=True, key="batch_wt")
    with batch_c2:
        batch_mut = st.file_uploader("Mutant 파일 여러 개", type=["gb", "gbk"],
                                     accept_multiple_files=True, key="batch_mut")
    if batch_wt and batch_mut:
        try:
            upload_lookup(batch_wt)
            upload_lookup(batch_mut)
        except ValueError as error:
            st.error(str(error))
        else:
            suggested = excel_config.copy() if excel_config is not None else suggest_batch_pairs(batch_wt, batch_mut, wt_label, mut_label)
            st.caption("파일명의 CNAG 번호 또는 공통 이름이 같은 짝을 추천합니다. WT 후보가 여러 개이거나 이름이 다르면 빈칸에서 직접 선택하세요. 같은 WT를 여러 Mutant와 짝지을 수 있습니다.")
            pair_editor_key = "batch_pairs_" + hashlib.sha256(
                json.dumps([[f.name for f in batch_wt], [f.name for f in batch_mut],
                            hashlib.sha256(excel_upload.getvalue()).hexdigest() if excel_upload else ""], ensure_ascii=False).encode("utf-8")
            ).hexdigest()[:12]
            batch_pairs = st.data_editor(suggested, hide_index=True, use_container_width=True,
                key=pair_editor_key, column_config={
                    "디자인 ID": st.column_config.TextColumn("디자인 ID", disabled=True),
                    "출력": st.column_config.CheckboxColumn("출력"),
                    "WT 파일": st.column_config.SelectboxColumn("WT 파일", options=list(dict.fromkeys([""] + [f.name for f in batch_wt] + suggested["WT 파일"].tolist())), required=True),
                    "Mutant 파일": st.column_config.SelectboxColumn("Mutant 파일",
                        options=list(dict.fromkeys([f.name for f in batch_mut] + suggested["Mutant 파일"].tolist())), required=True),
                    "출력 이름": st.column_config.TextColumn("출력 이름", required=True),
                    "WT 라벨": st.column_config.TextColumn("WT 라벨"),
                    "Mutant 라벨": st.column_config.TextColumn("Mutant 라벨"),
                })
            batch_formats = st.multiselect("출력 형식", ["PNG", "SVG", "PDF"],
                                           default=["PNG", "SVG", "PDF"], key="batch_formats")
            with st.expander("일괄 작업의 HR 점선 연결"):
                st.caption("아래 연결을 모든 파일 짝에 적용합니다. 해당 이름이 없는 그림에서는 그 연결을 건너뜁니다.")
                batch_hr = st.data_editor(default_hr_pairs(), num_rows="dynamic", hide_index=True,
                                          use_container_width=True, key="batch_hr")
            batch_settings = {key: st.session_state[key] for key in PRESET_DEFAULTS}
            signature = batch_signature(batch_wt, batch_mut, batch_pairs, batch_settings,
                                        batch_formats, batch_hr, rules_key, excel_styles)
            selected_pairs = batch_pairs[batch_pairs["출력"] == True]
            if st.button("📦 선택한 디자인 일괄 생성", key="batch_generate"):
                if selected_pairs.empty:
                    st.warning("출력할 항목을 하나 이상 선택해 주세요.")
                elif excel_upload and excel_config is None:
                    st.warning("엑셀 설정의 오류를 먼저 수정하거나 업로드를 해제해 주세요.")
                elif not batch_formats:
                    st.warning("출력 형식을 하나 이상 선택해 주세요.")
                elif any(str(name) not in {f.name for f in batch_wt} for name in selected_pairs["WT 파일"]):
                    st.warning("선택한 모든 행에서 업로드한 WT 파일의 짝을 지정해 주세요.")
                elif any(str(name) not in {f.name for f in batch_mut} for name in selected_pairs["Mutant 파일"]):
                    st.warning("선택한 모든 행의 Mutant 파일을 업로드하거나 파일명을 수정해 주세요.")
                else:
                    progress = st.progress(0.0, text="디자인 생성 중...")
                    try:
                        data, report, completed = build_batch_zip(
                            batch_pairs, batch_wt, batch_mut, rules_key, batch_settings,
                            batch_formats, batch_hr,
                            progress_callback=lambda done, total: progress.progress(
                                done / total, text=f"디자인 생성: {done}/{total}"), text_styles=excel_styles)
                        st.session_state["batch_result"] = (signature, data, report, completed)
                    except Exception as error:
                        st.error(f"일괄 작업 오류: {error}")
                    finally:
                        progress.empty()
            result = st.session_state.get("batch_result")
            if result and result[0] == signature:
                _, zip_data, report, completed = result
                if completed:
                    st.success(f"{completed}개 디자인 생성 완료. 각 디자인 폴더에 선택한 형식이 저장됩니다.")
                    st.download_button("📥 디자인 전체 ZIP 다운로드", data=zip_data,
                        file_name="allele_designs.zip", mime="application/zip", key="batch_download")
                if (report["상태"] == "실패").any():
                    st.warning("일부 파일을 처리하지 못했습니다. 아래 표의 내용을 확인해 주세요.")
                st.dataframe(report, hide_index=True, use_container_width=True)
    else:
        st.info("위의 WT·Mutant 업로드 영역에 각각 여러 GenBank 파일을 선택해 주세요.")


if not (wt_file and mutant_file):
    with tab_edit:
        st.info("먼저 '파일 업로드' 탭에서 GenBank 파일 두 개를 올려주세요.")
    with tab_hr:
        st.info("파일 업로드 후 사용할 수 있습니다.")
    with tab_view:
        st.info("파일 업로드 후 미리보기가 표시됩니다.")
    with tab_free:
        st.info("파일 업로드 후 사용할 수 있습니다.")
    st.stop()

single_styles = rules_for_design(excel_styles, "")
if excel_upload and excel_config is None:
    st.error("엑셀 설정의 오류를 수정하거나 엑셀 업로드를 해제해 주세요.")
    st.stop()
if excel_config is not None:
    matching = excel_config[(excel_config["WT 파일"] == wt_file.name) &
                            (excel_config["Mutant 파일"] == mutant_file.name)]
    if not matching.empty:
        with tab_upload:
            selected_design_id = st.selectbox("적용할 엑셀 디자인", matching["디자인 ID"].tolist(), key="single_excel_design")
            selected = matching[matching["디자인 ID"] == selected_design_id].iloc[0]
            st.caption("이 파일 짝의 엑셀 라벨과 단어별 서식을 적용합니다.")
        wt_label, mut_label = selected["WT 라벨"], selected["Mutant 라벨"]
        single_styles = rules_for_design(excel_styles, selected_design_id)
    else:
        with tab_upload:
            st.info("엑셀에 이 파일 짝이 없어서 기본 라벨과 공통 서식만 적용합니다.")

wt_df_raw, wt_len, wt_notes = parse_to_dataframe(
    decode_genbank_bytes(wt_file.getvalue()), False, rules_key)
mut_df_raw, mut_len, mut_notes = parse_to_dataframe(
    decode_genbank_bytes(mutant_file.getvalue()), True, rules_key)

editor_config = {
    "표시": st.column_config.CheckboxColumn("표시"),
    "이름": st.column_config.TextColumn("📝 이름"),
    "길이_배수": None,
    "Y_띄우기_mm": None,
    "위로향함": st.column_config.CheckboxColumn("⬆️ 위로"),
    "종류": st.column_config.SelectboxColumn("🏷️ 종류", options=CATEGORIES, required=True),
    "화살표": st.column_config.SelectboxColumn(
        "➡️ 화살표", options=ARROW_CHOICES, required=True,
        help="박스(Exon·Marker·Promoter)를 블록 화살표로 그립니다. 자동 = Marker는 각 strand 방향으로, Exon은 마지막 exon만 화살표로 표시합니다."),
    "시작": st.column_config.NumberColumn("시작(bp)", disabled=True),
    "종료": None, "방향": None,
}


def summarize(df, notes, title):
    st.markdown(f"**{title}**")
    counts = df["종류"].value_counts().to_dict() if not df.empty else {}
    shown = int(df["표시"].sum()) if not df.empty else 0
    st.caption(
        f"feature {len(df)}개 중 {shown}개 표시 · "
        + (" · ".join(f"{k} {v}" for k, v in counts.items()) if counts else "없음")
    )
    if notes:
        with st.expander(f"⚠️ 자동 추정한 항목 {len(notes)}개 — 확인 권장", expanded=False):
            for n in notes:
                st.write("•", n)

with tab_edit:
    st.caption(
        "GenBank의 모든 feature를 그대로 가져옵니다. 분류가 틀렸으면 **🏷️ 종류**를 드롭다운에서 "
        "바꾸면 그림에 즉시 반영됩니다. 여기서 정하는 건 네 가지뿐입니다 — "
        "**그릴지(표시) · 뭐라고 쓸지(이름) · 무슨 종류인지(종류) · 선의 위인지 아래인지(⬆️ 위로)**. "
        "세부 위치는 '✍️ 직접 편집' 탭에서 마우스로 옮기세요."
    )
    ec1, ec2 = st.columns(2)
    with ec1:
        summarize(wt_df_raw, wt_notes, "WT")
        wt_df = st.data_editor(wt_df_raw, column_config=editor_config, hide_index=True,
                               use_container_width=True, key="wt_ed")
    with ec2:
        summarize(mut_df_raw, mut_notes, "Mutant")
        mut_df = st.data_editor(mut_df_raw, column_config=editor_config, hide_index=True,
                                use_container_width=True, key="mut_ed")

with tab_hr:
    st.caption("편집한 '이름'을 입력하면 X자 점선으로 연결합니다.")
    default_hr = default_hr_pairs()
    hr_edited = st.data_editor(default_hr, num_rows="dynamic", hide_index=True, key="hr_editor")

figure_settings = {key: st.session_state[key] for key in PRESET_DEFAULTS}
fig, figure_info = render_gene_map(
    wt_df, mut_df, wt_len, mut_len, hr_edited, figure_settings, wt_label, mut_label, single_styles)
fig_height = figure_info["fig_height"]
dropped_wt, dropped_mut = figure_info["dropped_wt"], figure_info["dropped_mut"]
wt_long, mut_long = figure_info["wt_long"], figure_info["mut_long"]

# ---------- 미리보기 & 다운로드 ----------
with tab_view:
    st.caption(
        f"현재 출력 크기: **{fig_width:.1f} × {fig_height:.2f} inch** "
        f"({fig_width * 25.4:.0f} × {fig_height * 25.4:.0f} mm) · 세로는 요소 수에 맞춰 자동 계산됩니다."
    )
    if dropped_wt or dropped_mut:
        parts = []
        if dropped_wt:
            parts.append("WT: " + ", ".join(dropped_wt))
        if dropped_mut:
            parts.append("Mutant: " + ", ".join(dropped_mut))
        st.warning(
            "양 끝 자르기 범위 밖이라 그려지지 않은 요소가 있습니다 — "
            + " / ".join(parts)
            + " · 사이드바 '양 끝 자르기' 값을 줄이면 다시 나타납니다."
        )
    if wt_long or mut_long:
        st.info(
            "박스보다 이름이 길어 박스 밖에 표기한 항목: "
            + ", ".join(sorted(set(wt_long + mut_long)))
            + " · '요소 편집' 탭에서 이름을 줄이면 박스 안으로 들어갑니다."
        )
    preview = BytesIO()
    fig.savefig(preview, format="png", dpi=150, bbox_inches="tight",
                facecolor="white", pad_inches=0.05)
    st.image(preview.getvalue(), use_container_width=True)

    st.markdown("---")
    d1, d2, d3 = st.columns(3)
    for col, fmt, mime, label, fname in [
        (d1, "png", "image/png", "📸 PNG", "figure.png"),
        (d2, "svg", "image/svg+xml", "🎨 SVG (편집 가능)", "figure.svg"),
        (d3, "pdf", "application/pdf", "📄 PDF (벡터)", "figure.pdf"),
    ]:
        buf = BytesIO()
        kw = {"dpi": dpi_setting} if fmt == "png" else {}
        fig.savefig(buf, format=fmt, bbox_inches="tight", pad_inches=0.05, **kw)
        with col:
            st.download_button(label, data=buf.getvalue(), file_name=fname,
                               mime=mime, use_container_width=True)
    st.caption("SVG·PDF는 글자가 텍스트로 살아 있어 Illustrator에서 그대로 수정할 수 있습니다.")



# ---------- ✍️ 직접 편집: 브라우저에서 요소를 마우스로 옮긴다 ----------
with tab_free:
    st.caption(
        "PowerPoint처럼 요소를 끌어서 옮길 수 있습니다. 라벨·눈금·박스 아무거나 클릭한 뒤 드래그하거나, "
        "선택 후 방향키(Shift = 크게)로 밀면 됩니다. **설정을 바꾸면 그림이 다시 그려지면서 "
        "드래그한 위치는 초기화되니, 미세 조정은 맨 마지막에 하세요.**"
    )

    svg_buf = BytesIO()
    fig.savefig(svg_buf, format="svg", bbox_inches="tight", pad_inches=0.05)
    svg_src = svg_buf.getvalue().decode("utf-8")
    svg_src = svg_src[svg_src.index("<svg"):]

    stage_h = int(min(1000, 820 * fig_height / max(fig_width, 0.1))) + 110
    html = r"""
<style>
  body { margin:0; font-family: Arial, sans-serif; }
  #bar { display:flex; gap:6px; flex-wrap:wrap; align-items:center;
         padding:6px 2px 10px; font-size:13px; }
  #bar button { font:inherit; padding:5px 10px; border:1px solid #cfcfcf;
                background:#fff; border-radius:6px; cursor:pointer; }
  #bar button:hover { background:#f2f2f2; }
  #hint { color:#666; font-size:12px; margin-left:4px; }
  #stage { border:1px solid #e0e0e0; border-radius:8px; background:#fff;
           padding:8px; overflow:auto; }
  .sel { filter: drop-shadow(0 0 2.5px #2563eb); }
  svg [id] { cursor: default; }
</style>
<div id="bar">
  <button id="undo">되돌리기 (Ctrl+Z)</button>
  <button id="reset">전부 제자리로</button>
  <button id="dlsvg">SVG 저장</button>
  <button id="dlpng">PNG 저장</button>
  <button id="copy">SVG 복사</button>
  <span id="hint">요소를 클릭 → 드래그 또는 방향키</span>
</div>
<div id="stage">__SVG__</div>
<script>
(function(){
  const stage = document.getElementById('stage');
  const svg = stage.querySelector('svg');
  const vb = svg.viewBox.baseVal;
  svg.removeAttribute('width'); svg.removeAttribute('height');
  svg.style.width = '100%'; svg.style.height = 'auto';
  svg.style.touchAction = 'none';

  const items = Array.prototype.slice.call(svg.querySelectorAll('g[id]'))
    .filter(function(g){
      return /^(text_|line2d_|patch_|FancyArrow|PathCollection|PolyCollection)/.test(g.id)
             && g.parentNode && g.parentNode.id !== 'figure_1';
    });
  items.forEach(function(g){ g.style.cursor = 'move'; g.style.pointerEvents = 'all'; });

  function getT(el){
    const m = /translate\(\s*([-\d.]+)[ ,]\s*([-\d.]+)\s*\)/.exec(el.getAttribute('transform') || '');
    return m ? [parseFloat(m[1]), parseFloat(m[2])] : [0, 0];
  }
  function setT(el, x, y){ el.setAttribute('transform', 'translate(' + x + ',' + y + ')'); }
  function scale(){ return vb.width / svg.getBoundingClientRect().width; }

  let sel = null, drag = null;
  const undoStack = [];

  function select(el){
    if (sel) sel.classList.remove('sel');
    sel = el;
    if (sel) sel.classList.add('sel');
  }

  items.forEach(function(g){
    g.addEventListener('pointerdown', function(ev){
      ev.preventDefault(); ev.stopPropagation();
      select(g);
      const t = getT(g);
      drag = {el: g, x0: ev.clientX, y0: ev.clientY, tx: t[0], ty: t[1], moved: false};
      g.setPointerCapture(ev.pointerId);
    });
    g.addEventListener('pointermove', function(ev){
      if (!drag || drag.el !== g) return;
      const k = scale();
      const nx = drag.tx + (ev.clientX - drag.x0) * k;
      const ny = drag.ty + (ev.clientY - drag.y0) * k;
      if (!drag.moved){ undoStack.push({el: g, t: [drag.tx, drag.ty]}); drag.moved = true; }
      setT(g, nx, ny);
    });
    g.addEventListener('pointerup', function(ev){
      if (drag && drag.el === g) g.releasePointerCapture(ev.pointerId);
      drag = null;
    });
  });

  stage.addEventListener('pointerdown', function(){ select(null); });

  document.addEventListener('keydown', function(ev){
    if (ev.ctrlKey && ev.key.toLowerCase() === 'z'){ undo(); ev.preventDefault(); return; }
    if (!sel) return;
    const step = (ev.shiftKey ? 10 : 2) * (vb.width / 800);
    let dx = 0, dy = 0;
    if (ev.key === 'ArrowLeft') dx = -step;
    else if (ev.key === 'ArrowRight') dx = step;
    else if (ev.key === 'ArrowUp') dy = -step;
    else if (ev.key === 'ArrowDown') dy = step;
    else return;
    ev.preventDefault();
    const t = getT(sel);
    undoStack.push({el: sel, t: t});
    setT(sel, t[0] + dx, t[1] + dy);
  });

  function undo(){
    const last = undoStack.pop();
    if (last) setT(last.el, last.t[0], last.t[1]);
  }
  document.getElementById('undo').onclick = undo;
  document.getElementById('reset').onclick = function(){
    items.forEach(function(g){ g.removeAttribute('transform'); });
    undoStack.length = 0;
  };

  function serialize(){
    const clone = svg.cloneNode(true);
    clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    clone.setAttribute('width', vb.width + 'pt');
    clone.setAttribute('height', vb.height + 'pt');
    return new XMLSerializer().serializeToString(clone);
  }
  function download(blob, name){
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob); a.download = name;
    document.body.appendChild(a); a.click(); a.remove();
  }
  document.getElementById('dlsvg').onclick = function(){
    download(new Blob([serialize()], {type: 'image/svg+xml'}), 'figure_edited.svg');
  };
  document.getElementById('dlpng').onclick = function(){
    const url = URL.createObjectURL(new Blob([serialize()], {type: 'image/svg+xml;charset=utf-8'}));
    const img = new Image();
    img.onload = function(){
      const k = 4, c = document.createElement('canvas');
      c.width = vb.width * k; c.height = vb.height * k;
      const ctx = c.getContext('2d');
      ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, c.width, c.height);
      ctx.drawImage(img, 0, 0, c.width, c.height);
      URL.revokeObjectURL(url);
      c.toBlob(function(b){ download(b, 'figure_edited.png'); });
    };
    img.src = url;
  };
  document.getElementById('copy').onclick = function(){
    navigator.clipboard.writeText(serialize()).then(function(){
      document.getElementById('hint').textContent = 'SVG를 클립보드에 복사했습니다.';
    });
  };
})();
</script>
"""
    components.html(html.replace("__SVG__", svg_src), height=stage_h, scrolling=True)
    st.caption(
        "저장 버튼이 브라우저 정책으로 막히면 'SVG 복사'를 눌러 텍스트 편집기에 붙여넣고 "
        ".svg로 저장하세요. 더 자유롭게 다듬으려면 SVG를 Illustrator·Inkscape·PowerPoint에서 "
        "열고 그룹 해제하면 모든 도형과 글자가 그대로 편집됩니다."
    )


plt.close(fig)
# 세션이 계속 재실행돼도 Matplotlib Figure가 누적되지 않게 정리.
plt.close(fig)
