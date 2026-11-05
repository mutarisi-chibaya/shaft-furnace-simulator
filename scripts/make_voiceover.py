"""Voice-over for the demo video from docs/demo_script.md, with a South African English voice.

Usage:
    python scripts/make_voiceover.py --dry-run             # print the SSML, call nothing
    python scripts/make_voiceover.py [--voice en-ZA-LukeNeural] [--rate -5%]

Uses Microsoft Azure AI Speech (text to speech, REST API). Needs an Azure Speech resource;
set its key and region first (PowerShell):
    $env:AZURE_SPEECH_KEY = "<key>"
    $env:AZURE_SPEECH_REGION = "<region, e.g. southafricanorth or westeurope>"
The default voice en-ZA-LukeNeural is Azure's South African English male neural voice
(check the current voice list in the Azure Speech documentation if it has changed).

Narration = the lines starting with "> " under each "## Scene" heading. Writes
build/voiceover/scene_01.mp3 ... and demo_voiceover.mp3 (all scenes joined; MP3 frames of the
same format can be concatenated). The narration text is sent to Azure; nothing else is.
"""

from __future__ import annotations

import argparse
import os
import re
import urllib.request
from pathlib import Path
from xml.sax.saxutils import escape

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "docs" / "demo_script.md"
OUT = ROOT / "build" / "voiceover"
#: Pronunciation hints (written form -> spoken form) for words a TTS voice may misread.
ALIASES = {
    "wüstite": "woostite",
    "REDUCTOR": "re-duck-tor",
    "HyIron": "high-iron",
}


def scenes(path: Path) -> list[tuple[str, list[str]]]:
    """(scene title, narration paragraphs) in order of the script."""
    out: list[tuple[str, list[str]]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("## Scene"):
            out.append((line[3:].strip(), []))
        elif line.startswith("> ") and out:
            out[-1][1].append(line[2:].strip())
    return [s for s in out if s[1]]


def ssml(paragraphs: list[str], voice: str, rate: str) -> str:
    """SSML for one scene: paragraphs separated by short pauses, aliases applied."""
    parts = []
    for p in paragraphs:
        text = escape(p)
        for written, spoken in ALIASES.items():
            text = re.sub(
                rf"\b{re.escape(written)}\b", f'<sub alias="{spoken}">{written}</sub>', text
            )
        parts.append(text)
    body = '<break time="700ms"/>'.join(parts)
    return (
        '<speak version="1.0" xmlns="http://www.w3.org/2001/10/synthesis" xml:lang="en-ZA">'
        f'<voice name="{voice}"><prosody rate="{rate}">{body}</prosody></voice></speak>'
    )


def synthesize(text: str, key: str, region: str) -> bytes:
    req = urllib.request.Request(
        f"https://{region}.tts.speech.microsoft.com/cognitiveservices/v1",
        data=text.encode("utf-8"),
        headers={
            "Ocp-Apim-Subscription-Key": key,
            "Content-Type": "application/ssml+xml",
            "X-Microsoft-OutputFormat": "audio-24khz-96kbitrate-mono-mp3",
            "User-Agent": "shaftsim-voiceover",
        },
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=120) as r:
        return bytes(r.read())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--voice", default="en-ZA-LukeNeural")
    ap.add_argument("--rate", default="-5%", help="speaking rate, e.g. -5%% (slower) or +0%%")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    sc = scenes(SCRIPT)
    words = sum(len(p.split()) for _, ps in sc for p in ps)
    print(f"{len(sc)} scenes, {words} words (~{words / 150:.1f} min at 150 words/min)")
    if a.dry_run:
        for title, ps in sc:
            print(f"\n# {title}\n{ssml(ps, a.voice, a.rate)}")
        return
    key, region = os.environ.get("AZURE_SPEECH_KEY"), os.environ.get("AZURE_SPEECH_REGION")
    if not key or not region:
        raise SystemExit("set AZURE_SPEECH_KEY and AZURE_SPEECH_REGION (see the docstring)")
    OUT.mkdir(parents=True, exist_ok=True)
    joined = bytearray()
    for i, (title, ps) in enumerate(sc, start=1):
        audio = synthesize(ssml(ps, a.voice, a.rate), key, region)
        (OUT / f"scene_{i:02d}.mp3").write_bytes(audio)
        joined += audio
        print(f"  scene {i:02d} ({title}): {len(audio) / 1024:.0f} kB")
    (OUT / "demo_voiceover.mp3").write_bytes(bytes(joined))
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
