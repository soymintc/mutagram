import argparse
import os
import re
import glob
from PIL import Image
import matplotlib.pyplot as plt
import io
import igver
def resize_png(fig, out_path):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=100, bbox_inches='tight', pad_inches=0)
    buf.seek(0)
    img = Image.open(buf)
    cropped = img.crop(box=(100, 10, img.width-100, img.height-20))
    resized = cropped.resize((20, 50), Image.Resampling.LANCZOS)
    resized.save(out_path)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam_dir', required=True)
    parser.add_argument('--out_dir', required=True)
    parser.add_argument('--resized_dir', required=True)
    parser.add_argument('--region', required=True)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.resized_dir, exist_ok=True)

    for bam_file in glob.glob(os.path.join(args.bam_dir, '*.bam')):
        letter = os.path.basename(bam_file).split('.')[0]
        png_path = os.path.join(args.out_dir, f'{letter}.png')
        png_resized = os.path.join(args.resized_dir, f'{letter}.png')
        figs = igver.load_screenshots([bam_file], [args.region], max_panel_height=800, overlap_display='expand')
        fig = figs[0]
        fig.savefig(png_path)
        resize_png(fig, png_resized)

if __name__ == '__main__':
    main()
