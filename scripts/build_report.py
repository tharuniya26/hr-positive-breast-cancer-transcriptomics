"""Render the editable Quarto-flavoured Markdown report with ReportLab.

This lightweight renderer supports the document's headings, paragraphs and
standalone figure directives. It does not claim to implement Quarto itself.
"""
from pathlib import Path
import html,re
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib.units import cm
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Image,PageBreak,KeepTogether

source=Path("reports/molecular_transcriptomic_analysis.qmd")
target=Path("reports/molecular_transcriptomic_analysis.pdf")
raw=source.read_text(encoding="utf-8")
body=raw.split("---",2)[2].strip()
styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name="TitleCustom",parent=styles["Title"],fontName="Helvetica-Bold",fontSize=19,leading=24,spaceAfter=22,textColor=colors.HexColor("#173B56")))
styles.add(ParagraphStyle(name="SectionCustom",parent=styles["Heading1"],fontSize=13,leading=16,spaceBefore=16,spaceAfter=8,textColor=colors.HexColor("#173B56")))
styles.add(ParagraphStyle(name="BodyCustom",parent=styles["BodyText"],fontSize=9.5,leading=14,spaceAfter=9))
styles.add(ParagraphStyle(name="CaptionCustom",parent=styles["BodyText"],fontSize=8,leading=11,alignment=TA_CENTER,spaceAfter=14,textColor=colors.HexColor("#3D5262")))

def markup(text):
    text=html.escape(text).replace("&lt;", "&lt;")
    text=re.sub(r"\*([^*]+)\*",r"<i>\1</i>",text)
    text=re.sub(r"`([^`]+)`",r"<font name='Courier'>\1</font>",text)
    return text

story=[Spacer(1,4*cm),Paragraph("Molecular and Transcriptomic Analysis of Hormone Receptor-Positive Breast Cancer",styles["TitleCustom"]),
       Paragraph("Tharuniya Anpalakan",styles["Heading2"]),Paragraph("Computational research report · 16 October 2026",styles["BodyCustom"]),
       Spacer(1,2*cm),Paragraph("Secondary analysis of the public GSE45827 breast-cancer microarray cohort",styles["BodyCustom"]),PageBreak()]
for block in re.split(r"\n\s*\n",body):
    block=block.strip()
    if not block:continue
    if block.startswith("# "):
        story.append(Paragraph(markup(block[2:]),styles["SectionCustom"]))
    elif block.startswith("!["):
        m=re.fullmatch(r"!\[(.*?)\]\((.*?)\)",block)
        if not m:raise ValueError("Malformed figure directive")
        path=(source.parent/m.group(2)).resolve()
        if not path.is_file():raise FileNotFoundError(path)
        img=Image(str(path)); ratio=min((17*cm)/img.imageWidth,(11*cm)/img.imageHeight)
        img.drawWidth=img.imageWidth*ratio;img.drawHeight=img.imageHeight*ratio
        story.append(KeepTogether([img,Spacer(1,.2*cm),Paragraph(markup(m.group(1)),styles["CaptionCustom"])]))
    else:story.append(Paragraph(markup(block.replace("\n"," ")),styles["BodyCustom"]))

def decorate(canvas,doc):
    canvas.saveState();w,h=A4
    canvas.setStrokeColor(colors.HexColor("#91AABD"));canvas.line(2*cm,h-1.7*cm,w-2*cm,h-1.7*cm)
    canvas.setFont("Helvetica",8);canvas.setFillColor(colors.HexColor("#3D5262"))
    canvas.drawString(2*cm,1.25*cm,"Tharuniya Anpalakan · GSE45827 secondary analysis")
    canvas.drawRightString(w-2*cm,1.25*cm,str(doc.page))
    canvas.restoreState()

doc=SimpleDocTemplate(str(target),pagesize=A4,rightMargin=2*cm,leftMargin=2*cm,
                      topMargin=2.2*cm,bottomMargin=2*cm,title="Molecular and Transcriptomic Analysis of Hormone Receptor-Positive Breast Cancer",author="Tharuniya Anpalakan")
doc.build(story,onFirstPage=decorate,onLaterPages=decorate)
print(target,target.stat().st_size)
