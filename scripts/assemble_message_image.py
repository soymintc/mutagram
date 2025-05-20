import argparse
import os
import re
import glob
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

def load_char2path(png_dir):
    char2path = defaultdict(dict)
    for png_path in glob.glob(os.path.join(png_dir, '*.png')):
        fname = os.path.basename(png_path)
        match = re.match(r'(.+?)\.fg([ACGT])\.bg([ACGT])\.png', fname)
        if match:
            dst_letter, fg, bg = match.groups()
            key = fg + bg
            char2path[key][dst_letter] = png_path
    return char2path

def parse_message(message_path):
    if os.path.exists(message_path):
        with open(message_path) as f:
            return f.read().strip().split('\n')
    else:
        return message_path.split('\n')

def prepare_letter_map(chars):
    letter_map = {x: x for x in chars}
    letter_map['.'] = 'dot'
    letter_map[','] = 'comma'
    letter_map['('] = 'parenthesis0'
    letter_map[')'] = 'parenthesis1'
    letter_map['!'] = 'exclamation'
    return letter_map

def draw_message(lines, char2path, letter_map, red_key, blue_key, output, line_spacing=1.3):
    sample_img = mpimg.imread(next(iter(char2path[blue_key].values())))
    img_h, img_w = sample_img.shape[:2]
    n_rows = len(lines)
    n_cols = max(len(line) for line in lines)

    fig_w = img_w * n_cols / 70
    fig_h = img_h * n_rows / 70
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=100)
    ax.set_xlim(0, img_w * n_cols)
    ax.set_ylim(0, img_h * n_rows * line_spacing)
    ax.axis('off')

    write_red = False
    for row_idx, line in enumerate(lines[::-1]):
        col_idx = -1
        for ch in line:
            col_idx += 1
            if ch == ' ':
                continue
            if ch == '_':
                col_idx -= 1
                write_red = not write_red
                continue
            path_map = char2path[red_key] if write_red else char2path[blue_key]
            dst = letter_map.get(ch)
            path = path_map.get(dst)
            if not path:
                continue
            img = mpimg.imread(path)
            x0 = col_idx * img_w
            y0 = row_idx * img_h * line_spacing
            ax.imshow(img, extent=[x0, x0 + img_w, y0, y0 + img_h])

    fig.savefig(output, bbox_inches='tight')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--png_dir', required=True)
    parser.add_argument('--message', required=True)
    parser.add_argument('--fg', default='C')
    parser.add_argument('--bg', default='A')
    parser.add_argument('-o', '--output', required=True)
    args = parser.parse_args()

    blue_key = args.fg + args.bg
    red_key = 'TA'
    char2path = load_char2path(args.png_dir)
    lines = parse_message(args.message)
    letter_map = prepare_letter_map(set(''.join(lines)))
    draw_message(lines, char2path, letter_map, red_key, blue_key, args.output)

if __name__ == '__main__':
    main()
