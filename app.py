import streamlit as st
import pandas as pd
import re
from Bio import SeqIO, Restriction
from io import StringIO, BytesIO
import matplotlib.pyplot as plt
import matplotlib.patches as patches

st.set_page_config(layout="wide", page_title="논문 피규어 스튜디오")
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
# 🎛️ 사이드바 — 기능별 Expander 그룹화
# ============================================================
# --- 프리셋 기본값 정의 (저장/복원 대상) ---
import json
PRESET_DEFAULTS = {
    "label_offset_x": 150, "text_color": "#000000",
    "label_font_size": 14, "marker_font_size": 14,
    "enzyme_font_size": 14, "primer_font_size": 14,
    "color_exon": "#002060", "color_marker": "#FFC000",
    "color_enzyme": "#FF0000", "color_promoter": "#AAAAAA",
    "color_probe": "#92D050", "color_line": "#000000",
    "box_h": 0.05, "probe_h": 0.05, "enzyme_line_len": 0.05,
    "primer_offset_base": 0.07, "line_width": 2.5,
    "fig_width": 16, "fig_height": 5, "dpi_setting": 300,
    "zoom_x": 1.0,
    "trim_left": 0, "trim_right": 0,
    "show_scalebar": True, "scalebar_bp": 1000,
    "scalebar_color": "#000000", "scalebar_lw": 2.5,
    "show_hr": True,
}
# session_state 초기화 (최초 1회)
for k, v in PRESET_DEFAULTS.items():
    if k not in st.session_state:
        st.session_state[k] = v

with st.sidebar:
    st.header("⚙️ 컨트롤 패널")

    with st.expander("💾 프리셋 (내 기본값)", expanded=True):
        st.caption("현재 설정을 JSON으로 저장하거나, 저장한 파일을 불러와서 한 번에 복원합니다.")

        # Export
        current_preset = {k: st.session_state[k] for k in PRESET_DEFAULTS}
        st.download_button(
            "⬇️ 현재 설정을 JSON으로 저장",
            data=json.dumps(current_preset, indent=2, ensure_ascii=False),
            file_name="figure_preset.json",
            mime="application/json",
            use_container_width=True,
        )

        # Import
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
                    st.success("프리셋을 적용했습니다!")
                    st.rerun()
                except Exception as e:
                    st.error(f"프리셋 파일을 읽지 못했습니다: {e}")

        if st.button("🔄 공장 초기값으로 리셋", use_container_width=True):
            for k, v in PRESET_DEFAULTS.items():
                st.session_state[k] = v
            st.session_state.pop("_last_preset_sig", None)
            st.rerun()

    with st.expander("📝 라벨 & 텍스트", expanded=True):
        st.caption(
            "`**굵게**` · `*기울임*` · `^위첨자^` · `~아래첨자~`\n\n"
            "**조합 가능:** `***굵은기울임***`, `**car1^Δ^**`, `*p~TEF1~*`"
        )
        wt_label = st.text_input("WT 라벨", "", placeholder="예: WT")
        mut_label = st.text_input("Mutant 라벨", "", placeholder="예: *car1*Δ")
        label_offset_x = st.slider("라벨 X축 거리", 0, 1000, step=10, key="label_offset_x")
        text_color = st.color_picker("텍스트 색상", key="text_color")

    with st.expander("🔠 글씨 크기"):
        label_font_size = st.slider("라벨 (WT/Mut)", 6, 24, key="label_font_size")
        marker_font_size = st.slider("Marker", 6, 20, key="marker_font_size")
        enzyme_font_size = st.slider("Enzyme", 6, 20, key="enzyme_font_size")
        primer_font_size = st.slider("Primer", 6, 20, key="primer_font_size")

    with st.expander("🎨 색상"):
        color_exon = st.color_picker("Exon", key="color_exon")
        color_marker = st.color_picker("Marker", key="color_marker")
        color_enzyme = st.color_picker("Enzyme", key="color_enzyme")
        color_promoter = st.color_picker("Promoter", key="color_promoter")
        color_probe = st.color_picker("Probe", key="color_probe")
        color_line = st.color_picker("Backbone", key="color_line")

    with st.expander("📏 크기 & 두께"):
        box_h = st.slider("유전자 박스 두께", 0.01, 0.10, step=0.01, key="box_h")
        probe_h = st.slider("Probe 두께", 0.01, 0.10, step=0.01, key="probe_h")
        enzyme_line_len = st.slider("Enzyme 선 길이", 0.01, 0.10, step=0.01, key="enzyme_line_len")
        primer_offset_base = st.slider("Primer 선 길이", 0.05, 0.20, step=0.01, key="primer_offset_base")
        line_width = st.slider("Backbone 두께", 1.0, 5.0, step=0.5, key="line_width")

    with st.expander("🖼️ 캔버스"):
        st.caption(
            "**가로 확대 배율**을 키우면 글자·박스 두께는 그대로, "
            "유전자 영역만 옆으로 시원하게 펼쳐져 라벨 겹침이 풀립니다."
        )
        fig_width_base = st.slider("가로 (기본)", 5, 30, key="fig_width")
        zoom_x = st.slider("🔍 가로 확대 배율", 1.0, 10.0, step=0.5, key="zoom_x")
        fig_height = st.slider("세로", 3, 15, key="fig_height")
        dpi_setting = st.select_slider("DPI", options=[150, 300, 600], key="dpi_setting")
        # 실제 figure에 들어갈 가로 = 기본 × 배율
        fig_width = fig_width_base * zoom_x

    with st.expander("✂️ 양 끝 자르기 (WT/Mutant 동시 적용)"):
        st.caption(
            "주요 영역에 집중하기 위해 자동 계산된 범위에서 양쪽 끝을 추가로 잘라냅니다. "
            "0이면 비활성. 단위는 bp."
        )
        trim_left = st.number_input("⬅️ 왼쪽 잘라내기 (bp)", min_value=0, step=100, key="trim_left")
        trim_right = st.number_input("➡️ 오른쪽 잘라내기 (bp)", min_value=0, step=100, key="trim_right")

    with st.expander("📏 Scale bar"):
        show_scalebar = st.checkbox("Scale bar 표시", key="show_scalebar")
        scalebar_bp = st.number_input("길이 (bp)", min_value=100, step=100, key="scalebar_bp")
        scalebar_color = st.color_picker("색상", key="scalebar_color")
        scalebar_lw = st.slider("두께", 1.0, 6.0, step=0.5, key="scalebar_lw")
        st.caption("1000bp 이상은 자동으로 'kb' 단위로 표기됩니다.")

    with st.expander("🔧 기타"):
        show_hr = st.checkbox("HR 점선 표시", key="show_hr")

    st.markdown("---")
    st.caption(
        "⚠️ **As-is 제공** · 결과물 검증 책임은 사용자에게 있습니다.\n\n"
        "🔒 업로드 파일은 서버에 저장되지 않고 세션 종료 시 파기됩니다.\n\n"
        "📜 Noncommercial Academic Use / 비영리 학술 사용"
    )
    st.caption("👨‍💻 Developed by Yujin Lee | 💡 버그 제보 및 기능 건의 환영!")

# ============================================================
# 텍스트 파서 (버그 수정)
# ============================================================
def format_text(text):
    """
    중첩/조합 서식 지원:
      ***x*** → 굵은 기울임
      **x**   → 굵게
      *x*     → 기울임
      ^x^     → 위첨자
      ~x~     → 아래첨자
    조합 예: **car1^Δ^**, *p**TEF1*** 등
    """
    if not text:
        return ""
    text = str(text)
    if "$" in text:           # 사용자가 직접 LaTeX 입력한 경우 그대로
        return text

    original = text
    has_markup = bool(re.search(r'\*\*\*|\*\*|\*|\^.+?\^|~.+?~', text))
    if not has_markup:
        return original

    # 1) 안쪽부터 바깥쪽으로 재귀 변환 (가장 긴 토큰 우선)
    #    굵은 기울임 → 굵게 → 기울임 → 위/아래 첨자
    def convert(s):
        prev = None
        while prev != s:
            prev = s
            # ***bold italic***
            s = re.sub(r'\*\*\*([^*]+?)\*\*\*', r'\\mathbf{\\mathit{\1}}', s)
            # **bold**
            s = re.sub(r'\*\*([^*]+?)\*\*', r'\\mathbf{\1}', s)
            # *italic*  (앞뒤가 별표가 아닐 때만)
            s = re.sub(r'(?<!\*)\*([^*]+?)\*(?!\*)', r'\\mathit{\1}', s)
            # ^superscript^
            s = re.sub(r'\^([^\^]+?)\^', r'^{\1}', s)
            # ~subscript~
            s = re.sub(r'~([^~]+?)~', r'_{\1}', s)
        return s

    converted = convert(text)

    # 2) mathtext 안에서 일반 텍스트(영문/숫자/한글)도 보존되도록
    #    LaTeX 명령부 바깥의 공백만 \ 로 이스케이프
    #    간단히: 전체 공백 → \  (mathtext에서 안전)
    converted = converted.replace(" ", r"\ ")

    return f"${converted}$"

# ============================================================
# GenBank 파서
# ============================================================
@st.cache_data
def parse_to_dataframe(file_content, is_mutant=False):
    record = SeqIO.read(StringIO(file_content), "genbank")
    data = []
    try:
        ALL_ENZYMES = {str(e).upper() for e in Restriction.AllEnzymes}
    except Exception:
        ALL_ENZYMES = {'ECORV','PSTI','BAMHI','ECORI','HINDIII','XHOI','NOTI','PPUMI'}

    primer_keywords = {'SO','L1','L2','R1','R2','PL','PR','PO','SO2','ILP','IRP'}

    for f in record.features:
        name = ""
        for key in ['label','note','gene','name','ApEinfo_label']:
            if key in f.qualifiers:
                name = str(f.qualifiers[key][0]).strip()
                break
        if not name:
            continue

        upper = name.upper()
        no_space = upper.replace(" ", "")
        start, end = int(f.location.start), int(f.location.end)
        strand = f.location.strand if f.location.strand is not None else 1
        ftype = f.type.lower()
        category = "Other"

        if 'PROBE' in upper:
            name, category = "Probe", "Probe"
        elif ftype == 'restriction_site' or upper in ALL_ENZYMES or no_space in ALL_ENZYMES:
            category = "Enzyme"
        elif any(m in upper for m in ['NAT','NEO','HYG','TRP1']):
            category = "Marker"
        elif 'PROMOTER' in upper:
            category = "Promoter"
        elif ftype in ['primer_bind','misc_binding'] or upper in primer_keywords or 'PRIMER' in upper:
            category = "Primer"
        elif ftype in ['cds','exon'] or upper.startswith('E'):
            category = "Exon"

        if category != "Other":
            # Mutant 파일의 Enzyme은 기본적으로 위로 향하도록
            up_default = (is_mutant and category == "Enzyme")
            data.append({
                "표시": upper not in ['QLP','QRP'],
                "이름": name, "길이_배수": 1.0, "Y축_띄우기": 0.0,
                "위로향함": up_default, "종류": category,
                "시작": start, "종료": end, "방향": strand
            })
    return pd.DataFrame(data), len(record)

# ============================================================
# 메인 화면 — 탭 구조
# ============================================================
tab_upload, tab_edit, tab_hr, tab_view = st.tabs(
    ["📂 파일 업로드", "🛠️ 요소 편집", "🔗 HR 연결", "🖼️ 미리보기 & 다운로드"]
)

with tab_upload:
    c1, c2 = st.columns(2)
    with c1:
        wt_file = st.file_uploader("WT GenBank (.gb)", type=["gb","gbk"])
    with c2:
        mutant_file = st.file_uploader("Mutant GenBank (.gb)", type=["gb","gbk"])

if wt_file and mutant_file:
    wt_df_raw, wt_len = parse_to_dataframe(wt_file.getvalue().decode("utf-8"), is_mutant=False)
    mut_df_raw, mut_len = parse_to_dataframe(mutant_file.getvalue().decode("utf-8"), is_mutant=True)

    editor_config = {
        "표시": st.column_config.CheckboxColumn("표시"),
        "이름": st.column_config.TextColumn("📝 이름"),
        "길이_배수": st.column_config.NumberColumn("📏 선 길이", min_value=0.1, max_value=5.0, step=0.1, format="%.1f"),
        "Y축_띄우기": st.column_config.NumberColumn("↕️ Y 띄우기", min_value=-0.5, max_value=0.5, step=0.01, format="%.2f"),
        "위로향함": st.column_config.CheckboxColumn("⬆️ 위로"),
        "종류": st.column_config.TextColumn("종류", disabled=True),
        "시작": None, "종료": None, "방향": None
    }

    with tab_edit:
        ec1, ec2 = st.columns(2)
        with ec1:
            st.markdown("**WT**")
            wt_df_edited = st.data_editor(wt_df_raw, column_config=editor_config, hide_index=True, key="wt_ed")
        with ec2:
            st.markdown("**Mutant**")
            mut_df_edited = st.data_editor(mut_df_raw, column_config=editor_config, hide_index=True, key="mut_ed")

    with tab_hr:
        st.caption("편집한 '이름'을 입력하면 X자 점선으로 연결합니다.")
        default_hr = pd.DataFrame({"WT_이름":["L1","L2","R1","R2"], "Mutant_이름":["L2","L1","R2","R1"]})
        hr_edited = st.data_editor(default_hr, num_rows="dynamic", hide_index=True, key="hr_editor")

    # ---------- 그리기 ----------
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Roboto','Arial','sans-serif']
    plt.rcParams['mathtext.default'] = 'regular'

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    y_wt, y_mut = 0.65, 0.35
    z_line, z_hr, z_primer, z_box = 1, 2, 4, 10

    # 표시되는 요소들의 실제 좌표 범위 자동 계산 → 빈 backbone 제거
    all_visible = pd.concat([
        wt_df_edited[wt_df_edited['표시']==True],
        mut_df_edited[mut_df_edited['표시']==True]
    ])
    if not all_visible.empty:
        data_start = int(all_visible['시작'].min())
        data_end   = int(all_visible['종료'].max())
    else:
        data_start, data_end = 0, max(wt_len, mut_len)

    # 좌우에 데이터 영역의 2%만큼 작은 패딩 추가 (최소 50bp)
    pad = max(int((data_end - data_start) * 0.02), 50)
    auto_start = max(0, data_start - pad)
    auto_end   = data_end + pad

    # 사용자 trim 적용 (양 끝에서 추가로 잘라냄)
    final_start = auto_start + trim_left
    final_end   = auto_end - trim_right
    if final_end <= final_start:   # 너무 많이 자른 경우 안전장치
        final_end = final_start + 100
    span = final_end - final_start

    def X(x):
        return x

    ax.plot([final_start, final_end],[y_wt,y_wt], color=color_line, lw=line_width, zorder=z_line)
    ax.plot([final_start, final_end],[y_mut,y_mut], color=color_line, lw=line_width, zorder=z_line)
    ax.text(final_start - label_offset_x, y_wt, format_text(wt_label), ha='right', va='center', color=text_color, fontsize=label_font_size)
    ax.text(final_start - label_offset_x, y_mut, format_text(mut_label), ha='right', va='center', color=text_color, fontsize=label_font_size)

    def draw(y_pos, df, is_top):
        for idx, row in df[df['표시']==True].iterrows():
            s,e,cat = row['시작'], row['종료'], row['종류']
            # 잘린 영역 밖이면 완전히 스킵
            if e <= final_start or s >= final_end: continue
            # Enzyme/Primer는 시작점이 영역 안에 있을 때만
            if cat in ("Enzyme", "Primer"):
                if s < final_start or s > final_end: continue
            # 박스/Probe는 영역 경계에 클램프
            s_clamp = max(s, final_start)
            e_clamp = min(e, final_end)
            sx, ex = s_clamp, e_clamp
            w = ex - sx
            sx_point = s
            fname = format_text(row['이름'])
            level, gap = row['길이_배수'], row['Y축_띄우기']
            d = 1 if row['위로향함'] else -1

            if cat == "Enzyme":
                y0 = y_pos + gap; y1 = y0 + enzyme_line_len*level*d
                ax.plot([sx_point,sx_point],[y0,y1], color=color_enzyme, lw=2, zorder=3)
                ax.text(sx_point, y1+0.02*d, fname, color=color_enzyme, ha='center',
                        va='bottom' if d==1 else 'top', fontsize=enzyme_font_size)
            elif cat == "Exon":
                ax.add_patch(patches.Rectangle((sx,y_pos-box_h/2), w, box_h, facecolor=color_exon, zorder=z_box))
            elif cat == "Marker":
                ax.add_patch(patches.Rectangle((sx,y_pos-box_h/2), w, box_h, facecolor=color_marker, zorder=z_box))
                ax.text(sx+w/2, y_pos, fname, ha='center', va='center', fontsize=marker_font_size, zorder=z_box+1)
            elif cat == "Promoter":
                ax.add_patch(patches.Rectangle((sx,y_pos-box_h/2), w, box_h, facecolor=color_promoter, zorder=z_box))
                ax.text(sx+w/2, y_pos, fname, ha='center', va='center', fontsize=marker_font_size, color='#000000', zorder=z_box+1)
            elif cat == "Probe":
                yp = y_pos+0.15 if y_pos==y_wt else y_pos-0.18
                ax.add_patch(patches.Rectangle((sx,yp), w, probe_h, facecolor=color_probe, alpha=0.8))
                ax.text(sx+w/2, yp+probe_h/2, fname, ha='center', va='center', fontsize=label_font_size)
            elif cat == "Primer":
                p = 1 if is_top else -1
                if row['위로향함']: p = -p
                y0 = y_pos + gap*p; y1 = y0 + primer_offset_base*p*level
                ax.plot([sx_point,sx_point],[y0,y1], color=color_line, lw=1, zorder=z_primer)
                dx = 120 if row['방향']==1 else -120
                ax.arrow(sx_point, y1, dx, 0, head_width=0.012, head_length=50,
                         fc=color_line, ec=color_line, lw=1, zorder=z_primer)
                ax.text(sx_point+dx/2, y1+0.015*p, fname, ha='center',
                        va='bottom' if p==1 else 'top', color=text_color, fontsize=primer_font_size, zorder=z_primer)

    draw(y_wt, wt_df_edited, True)
    draw(y_mut, mut_df_edited, False)

    if show_hr:
        for _, row in hr_edited.iterrows():
            wt_k = str(row.get('WT_이름','')).strip()
            mut_k = str(row.get('Mutant_이름','')).strip()
            if not wt_k or not mut_k or wt_k.lower()=='nan' or mut_k.lower()=='nan': continue
            wm = wt_df_edited[(wt_df_edited['표시']==True) & (wt_df_edited['이름']==wt_k)]
            mm = mut_df_edited[(mut_df_edited['표시']==True) & (mut_df_edited['이름']==mut_k)]
            if not wm.empty and not mm.empty:
                wx, mx = wm.iloc[0]['시작'], mm.iloc[0]['시작']
                if final_start <= wx <= final_end and final_start <= mx <= final_end:
                    ax.plot([X(wx), X(mx)],[y_wt,y_mut],'--', color=color_line, lw=1, alpha=0.5, zorder=z_hr)

    # Scale bar (왼쪽 위, WT 트랙 위)
    if show_scalebar:
        sb_x_start = final_start
        sb_x_end = final_start + scalebar_bp
        sb_y = 1.00   # WT(0.65)에서 충분히 위로 떨어뜨림
        ax.plot([sb_x_start, sb_x_end], [sb_y, sb_y],
                color=scalebar_color, lw=scalebar_lw, solid_capstyle='butt', zorder=20)
        # 단위 자동 변환
        if scalebar_bp >= 1000:
            val = scalebar_bp / 1000
            sb_label = f"{val:g} kb"   # :g → 1.0 → 1, 2.5 → 2.5
        else:
            sb_label = f"{scalebar_bp} bp"
        ax.text((sb_x_start + sb_x_end) / 2, sb_y + 0.02, sb_label,
                ha='center', va='bottom', color=scalebar_color,
                fontsize=label_font_size, zorder=20)

    # xlim: 자동 계산된 데이터 영역에 딱 맞춤. 좌측은 라벨용 여백만.
    ax.set_xlim(final_start - label_offset_x * 1.2, final_end)
    ax.set_ylim(0.1, 1.10)   # scale bar 공간 확보 (위쪽 여유 있게)
    ax.axis('off')

    with tab_view:
        import base64
        preview_buf = BytesIO()
        fig.savefig(preview_buf, format="png", dpi=100)
        preview_b64 = base64.b64encode(preview_buf.getvalue()).decode()

        col_mode, col_zoom = st.columns([1, 3])
        with col_mode:
            view_mode = st.radio("보기 모드", ["화면 맞춤", "확대"], horizontal=True, label_visibility="collapsed")
        with col_zoom:
            if view_mode == "확대":
                preview_zoom = st.slider("확대 (%)", 50, 400, 150, step=10, label_visibility="collapsed")
            else:
                preview_zoom = None

        if view_mode == "화면 맞춤":
            img_style = "display:block; max-width:100%; height:auto; margin:0 auto;"
            box_style = "overflow:hidden;"
        else:
            img_style = f"display:block; max-width:none; width:{preview_zoom}%; height:auto;"
            box_style = "overflow-x:auto; overflow-y:auto; max-height:600px;"

        st.markdown(
            f"""
            <div style="{box_style} border:1px solid #ddd; border-radius:6px; padding:8px; background:#fff;">
                <img src="data:image/png;base64,{preview_b64}" style="{img_style}" />
            </div>
            """,
            unsafe_allow_html=True,
        )
        st.caption("💡 '화면 맞춤'은 자동 축소, '확대'는 슬라이더로 배율을 조절하고 스크롤로 탐색합니다. 다운로드는 항상 풀 해상도.")
        st.markdown("---")
        d1,d2,d3 = st.columns(3)
        for col, fmt, mime, label, fname in [
            (d1,"png","image/png","📸 PNG","figure.png"),
            (d2,"svg","image/svg+xml","🎨 SVG","figure.svg"),
            (d3,"pdf","application/pdf","📄 PDF","figure.pdf"),
        ]:
            buf = BytesIO()
            kw = {"dpi":dpi_setting} if fmt=="png" else {}
            fig.savefig(buf, format=fmt, bbox_inches='tight', **kw)
            with col:
                st.download_button(label, data=buf.getvalue(), file_name=fname, mime=mime)
else:
    with tab_edit: st.info("먼저 '파일 업로드' 탭에서 GenBank 파일 두 개를 올려주세요.")
    with tab_hr: st.info("파일 업로드 후 사용할 수 있습니다.")
    with tab_view: st.info("파일 업로드 후 미리보기가 표시됩니다.")
