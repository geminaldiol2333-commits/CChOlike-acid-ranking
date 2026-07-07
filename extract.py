import docx, glob
doc = docx.Document(glob.glob(r'C:\Users\ycc\.agents\skills\ccho-acid-ranking\~\Downloads\CChO*.docx')[0])
with open('text.txt', 'w', encoding='utf-8') as f:
    f.write('\n'.join([p.text for p in doc.paragraphs]))
