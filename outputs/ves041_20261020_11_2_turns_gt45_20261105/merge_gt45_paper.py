"""Append the independently sized second paper page to the first page."""
import sys
from pathlib import Path
from pypdf import PdfReader, PdfWriter

first, second = map(Path, sys.argv[1:3])
writer = PdfWriter()
for source in (first, second):
    for page in PdfReader(source).pages:
        writer.add_page(page)
temporary = first.with_suffix('.merged.pdf')
with temporary.open('wb') as output:
    writer.write(output)
temporary.replace(first)
