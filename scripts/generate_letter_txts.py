import argparse
import os
import numpy as np
from PIL import Image, ImageDraw, ImageFont

def create_letter_image(letter, img_size=(120, 250), font_size=110, ttf="RobotoMono-Bold.ttf"):
    img = Image.new('L', img_size, color=255)
    draw = ImageDraw.Draw(img)
    font = ImageFont.truetype(ttf, font_size)
    text_bbox = draw.textbbox((0, 0), letter, font=font)
    text_w, text_h = text_bbox[2] - text_bbox[0], text_bbox[3] - text_bbox[1]
    position = ((img_size[0] - text_w) // 2, (img_size[1] - text_h) // 2)
    draw.text(position, letter, fill=0, font=font)
    return img

def resize_image(img, new_size):
    return img.resize(new_size, resample=Image.Resampling.BILINEAR)

def image_to_char_matrix(img, fg_char='C', bg_char='A', threshold=128, trim_top=15, trim_side=1):
    arr = np.array(img)
    matrix = np.where(arr > threshold, bg_char, fg_char)
    return matrix[trim_top:, trim_side:-trim_side]

def collapse_char_matrix(matrix):
    return '\n'.join(''.join(row) for row in matrix)

def create_letter_files(letters, letter_map, out_dir, fg='C', bg='A', resize=(70, 60), trim_top=20, trim_side=10, ttf="RobotoMono-Bold.ttf"):
    os.makedirs(out_dir, exist_ok=True)
    for letter in letters:
        dst_letter = letter_map.get(letter, letter)
        base_img = create_letter_image(letter, ttf=ttf)
        resized_img = resize_image(base_img, resize)
        char_matrix = image_to_char_matrix(resized_img, fg_char=fg, bg_char=bg, trim_top=trim_top, trim_side=trim_side)
        character = collapse_char_matrix(char_matrix)
        with open(f'{out_dir}/{dst_letter}.fg{fg}.bg{bg}.txt', 'w') as out:
            out.write(character + '\n')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out_dir', required=True)
    parser.add_argument('--font', default="RobotoMono-Bold.ttf")
    parser.add_argument('--fg', default='C')
    parser.add_argument('--bg', default='A')
    parser.add_argument('--letters', default='()')
    args = parser.parse_args()

    letters = list(args.letters)
    letter_map = {x: x for x in letters}
    letter_map['.'] = 'dot'
    letter_map[','] = 'comma'
    letter_map['('] = 'parenthesis0'
    letter_map[')'] = 'parenthesis1'
    letter_map['!'] = 'exclamation'

    create_letter_files(letters, letter_map, args.out_dir, fg=args.fg, bg=args.bg, ttf=args.font)

if __name__ == '__main__':
    main()
