"""Build an animated GIF of the Mirror-Score workflow from the static figure.

Frames: five progressive stage reveals (stages 1..k full opacity, remaining
stages dimmed) plus a final full-color frame. Assumes the five workflow
stages are laid out left-to-right in the source figure.
"""
import sys

from PIL import Image, ImageEnhance

SRC = "/workspace/mirrorscore/figures/fig_workflow.png"
OUT = "/workspace/mirrorscore/figures/fig_workflow.gif"
N_STAGES = 5
FRAME_MS = 1200
TARGET_W = 900


def build():
    img = Image.open(SRC).convert("RGB")
    w, h = img.size
    scale = TARGET_W / w
    img = img.resize((TARGET_W, int(h * scale)), Image.LANCZOS)
    w, h = img.size

    dim = ImageEnhance.Brightness(img).enhance(0.25)
    frames = []
    for k in range(1, N_STAGES + 1):
        x_split = int(w * k / N_STAGES)
        frame = Image.new("RGB", (w, h))
        frame.paste(dim, (0, 0))
        revealed = img.crop((0, 0, x_split, h))
        frame.paste(revealed, (0, 0))
        frames.append(frame)

    # last reveal frame is already the full figure; hold it longer
    durations = [FRAME_MS] * (len(frames) - 1) + [2500]
    frames[0].save(
        OUT, save_all=True, append_images=frames[1:], duration=durations,
        loop=0, optimize=True,
    )
    print(f"saved {OUT}: {len(frames)} frames, {w}x{h}")


if __name__ == "__main__":
    sys.exit(build())
