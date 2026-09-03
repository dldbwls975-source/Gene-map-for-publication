import json
import re
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

for _k, _v in PRESET_DEFAULTS.items():
    if _k not in st.session_state:
        st.session_state[_k] = _v


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
        try:
            with open("🧬 논문 피규어용 유전자 맵 스튜디오 사용 가이드.docx", "rb") as guide_file:
                st.download_button(
                    label="📄 가이드 다운로드 (.docx)",
                    data=guide_file,
                    file_name="논문_피규어_유전자맵_스튜디오_사용가이드.docx",
                    mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
                    use_container_width=True,
                )
        except FileNotFoundError:
            st.caption("가이드 파일이 서버에 없습니다.")

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
        wt_label = st.text_input("WT 라벨", "", placeholder="예: WT")
        mut_label = st.text_input("Mutant 라벨", "", placeholder="예: *car1*Δ")
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
def format_text(text):
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
            "길이_배수": 1.0,
            "Y_띄우기_mm": 0.0,
            "위로향함": up_default,
            "시작": start, "종료": end, "방향": strand,
        })

    df = pd.DataFrame(rows)
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
            y0 = base + 1.2 * MM
            y1 = y0 + cfg["probe_h"]
            any_outside = False
            for i, r in rows:
                w_bp = float(r["종료"]) - float(r["시작"])
                need = text_bp_width(r["이름"], fs, cfg["bp_per_inch"])
                lab_out = need > w_bp * 0.92
                any_outside = any_outside or lab_out
                plan[i] = {"kind": "Probe", "side": side, "y0": y0, "y1": y1,
                           "label": True, "fontsize": fs, "label_outside": lab_out}
            base = y1 + 1.0 * MM + (th + 0.8 * MM if any_outside else 0.0)

        extent[side] = max(extent[side], base)

    # --- 박스류: 라벨이 박스보다 길면 6pt까지 축소, 그래도 넘치면 스택 바깥으로 ---
    overflow, long_names = [], []
    for i, r in vis.iterrows():
        if r["종류"] not in ("Exon", "Marker", "Promoter"):
            continue
        w = float(r["종료"]) - float(r["시작"])
        fs = cfg["marker_font_size"]
        outside = False
        if r["종류"] in ("Marker", "Promoter"):
            need = text_bp_width(r["이름"], fs, cfg["bp_per_inch"])
            if need > w * 0.92:
                fs = max(6.0, fs * (w * 0.92) / max(need, 1e-9))
                if text_bp_width(r["이름"], fs, cfg["bp_per_inch"]) > w * 0.92:
                    outside = True
                    fs = cfg["marker_font_size"]
                    overflow.append((i, float(r["시작"]) + w / 2, str(r["이름"])))
                    long_names.append(str(r["이름"]))
        plan[i] = {"kind": r["종류"], "side": cfg["out_side"], "fontsize": fs,
                   "outside": outside}

    if overflow:
        side_out = 1 if extent[1] <= extent[-1] else -1
        th = text_height(cfg["marker_font_size"])
        y_out = extent[side_out] + 1.2 * MM
        placed = place_labels(overflow, cfg["marker_font_size"], cfg["bp_per_inch"],
                              gap_bp, max_shift, True, max_tiers=2)
        for i, x, name in overflow:
            tier, xlab = placed.get(i, (0, x))
            plan[i].update(side=side_out, y_out=y_out + tier * (th + 0.6 * MM), xlab=xlab)
        extent[side_out] = y_out + (max(t for t, _ in placed.values()) if placed else 0) \
            * (th + 0.6 * MM) + th

    return plan, extent, long_names


# ============================================================
# 메인
# ============================================================
tab_upload, tab_edit, tab_hr, tab_view, tab_free = st.tabs(
    ["📂 파일 업로드", "🛠️ 요소 편집", "🔗 HR 연결", "🖼️ 미리보기 & 다운로드", "✍️ 직접 편집"]
)

with tab_upload:
    c1, c2 = st.columns(2)
    with c1:
        wt_file = st.file_uploader("WT GenBank (.gb)", type=["gb", "gbk"])
    with c2:
        mutant_file = st.file_uploader("Mutant GenBank (.gb)", type=["gb", "gbk"])

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

rules_key = (
    split_rule(rule_marker), split_rule(rule_promoter),
    split_rule(rule_probe), split_rule(rule_primer),
    int(rule_enzyme_max), int(rule_primer_max), bool(hide_unnamed),
)
wt_df_raw, wt_len, wt_notes = parse_to_dataframe(
    wt_file.getvalue().decode("utf-8"), False, rules_key)
mut_df_raw, mut_len, mut_notes = parse_to_dataframe(
    mutant_file.getvalue().decode("utf-8"), True, rules_key)

editor_config = {
    "표시": st.column_config.CheckboxColumn("표시"),
    "이름": st.column_config.TextColumn("📝 이름"),
    "길이_배수": None,
    "Y_띄우기_mm": None,
    "위로향함": st.column_config.CheckboxColumn("⬆️ 위로"),
    "종류": st.column_config.SelectboxColumn("🏷️ 종류", options=CATEGORIES, required=True),
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
    default_hr = pd.DataFrame({"WT_이름": ["L1", "L2", "R1", "R2"],
                               "Mutant_이름": ["L2", "L1", "R2", "R1"]})
    hr_edited = st.data_editor(default_hr, num_rows="dynamic", hide_index=True, key="hr_editor")

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
    dedup_enzyme=dedup_enzyme, auto_stagger=auto_stagger,
)

wt_plan, wt_ext, wt_long = plan_track(wt_df, dict(cfg_common, out_side=1))
mut_plan, mut_ext, mut_long = plan_track(mut_df, dict(cfg_common, out_side=-1))

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

ax.text(final_start - 2 * MM * bp_per_inch, y_wt, format_text(wt_label),
        ha="right", va="center", color=text_color, fontsize=label_font_size)
ax.text(final_start - 2 * MM * bp_per_inch, y_mut, format_text(mut_label),
        ha="right", va="center", color=text_color, fontsize=label_font_size)


def draw_track(track_y, df, plan):
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
        name = format_text(row["이름"])

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
            if geo.get("label_outside"):
                y_lab = (hi if side == 1 else lo) + side * 0.8 * MM
                ax.text(sx + w / 2, y_lab, name, ha="center",
                        va="bottom" if side == 1 else "top",
                        fontsize=fs, color=text_color, zorder=Z_MARK)
            else:
                ax.text(sx + w / 2, (lo + hi) / 2, name, ha="center", va="center",
                        fontsize=fs, color=text_color, zorder=Z_MARK)

        elif kind in ("Exon", "Marker", "Promoter"):
            face = {"Exon": color_exon, "Marker": color_marker, "Promoter": color_promoter}[kind]
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
                ax.text(sx + w / 2, track_y, name, ha="center", va="center",
                        fontsize=fs, color=inner, zorder=Z_MARK)
    return dropped


dropped_wt = draw_track(y_wt, wt_df, wt_plan)
dropped_mut = draw_track(y_mut, mut_df, mut_plan)

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
