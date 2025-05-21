s a playful genomic-art pipeline that encodes characters as DNA sequences, simulates their alignment in `.bam` files, and visualizes them using IGV-like screenshots. The final output is a rendered image message composed entirely of genomic-style alignments.

![example](docs/example_output.png)

## 📦 Overview

The project consists of 4 modular Python scripts:

1. **`generate_letter_txts.py`** – Creates text representations of letters (`A`/`C` matrix) from a monospaced font.
2. **`make_letter_bams.py`** – Encodes these matrices as DNA reads and writes them into `.bam` files.
3. **`generate_letter_pngs.py`** – Uses `igver` (IGV automation) to generate screenshot PNGs from the BAMs and resizes them.
4. **`assemble_message_image.py`** – Composes resized PNGs into a message and displays it inline using matplotlib.

---

## 🚀 Quickstart

### 1. Generate letter `.txt` files

```bash
python generate_letter_txts.py \
  --letters '()!' \
  --out_dir results/letters \
  --font data/fonts/RobotoMono-Bold.ttf \
  --fg T --bg A
```

### 2. Create `.bam` files

```bash
python make_letter_bams.py \
  --letters_dir results/letters \
  --out_dir results/bam \
  --fasta path/to/GRCh37-lite.fa \
  --region_yaml metadata/regions.yaml \
  --fg T --bg A
```

### 3. Generate and resize PNGs

```bash
python generate_letter_pngs.py \
  --bam_dir results/bam \
  --out_dir results/png \
  --resized_dir results/png_resize \
  --region '1:100000-100070'
```

### 4. Render a message

At this point `--png_dir` should contain resized letter png files.

```bash
python assemble_message_image.py \
  --png_dir results/png_resize \
  --message messages/Foobar.txt \
  --fg T --bg A
```

Or

```bash
python assemble_message_image.py \
  --png_dir results/png_resize \
  --message "Hello World!" \
  --fg T --bg A
```

> Use underscores `_like this_` to toggle alternate colors (blue ↔ red) for words.

---

## 🧠 Requirements

- Python 3.7+
- Packages:
  - `pysam`, `pyfaidx`, `Pillow`, `matplotlib`, `yaml`, `igver`
- External:
  - `samtools` (must be in `$PATH`)

---

## 📁 Project Structure

```
mutagram/
├── scripts/
│   ├── generate_letter_txts.py
│   ├── make_letter_bams.py
│   ├── generate_letter_pngs.py
│   └── assemble_message_image.py
├── metadata/
│   └── regions.yaml
├── messages/
│   └── mymessage.txt
├── results/
│   ├── letters/
│   ├── bam/
│   ├── png/
│   └── png_resize/
├── data/
│   └── fonts/
│       └── RobotoMono-Bold.ttf
```

---

## 📝 License

This project is licensed under the [MIT License](LICENSE).

---

## 👨‍🔬 Author

Created by Seongmin Choi @ Memorial Sloan Kettering.

If you enjoy this or want to build creative genomic artwork, feel free to contribute or fork the repo!

a
