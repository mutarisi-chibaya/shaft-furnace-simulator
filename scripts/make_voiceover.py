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
import shutil
import subprocess
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


def plain_text(paragraphs: list[str]) -> str:
    """Narration as plain text with the pronunciation aliases applied (offline engine)."""
    text = " ".join(paragraphs)
    for written, spoken in ALIASES.items():
        text = re.sub(rf"{re.escape(written)}", spoken, text)
    return text


def synthesize_windows(text: str, voice: str, wav: Path) -> None:
    """Offline speech with an installed Windows voice (System.Speech via PowerShell)."""
    ps = (
        "Add-Type -AssemblyName System.Speech;"
        "$s = New-Object System.Speech.Synthesis.SpeechSynthesizer;"
        f"$s.SelectVoice('{voice}'); $s.Rate = -1;"
        f"$s.SetOutputToWaveFile('{wav}');"
        "$s.Speak([Console]::In.ReadToEnd()); $s.Dispose()"
    )
    subprocess.run(["powershell", "-NoProfile", "-Command", ps], input=text, text=True, check=True)


def synthesize_edge(text: str, voice: str, rate: str, mp3: Path) -> None:
    """Neural voice through the edge-tts package (Microsoft Edge's online read-aloud voices).

    Unofficial: edge-tts imitates the Edge browser's read-aloud client. It needs no key but
    may stop working if Microsoft changes the service, and it is not licensed for commercial
    use. Plain text only (no SSML), so the pronunciation aliases are applied as text.
    """
    import asyncio

    import edge_tts

    asyncio.run(edge_tts.Communicate(text, voice=voice, rate=rate).save(str(mp3)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--voice", default="en-ZA-LukeNeural")
    ap.add_argument("--rate", default="-5%", help="speaking rate, e.g. -5%% (slower) or +0%%")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument(
        "--engine",
        choices=("azure", "edge", "windows"),
        default="azure",
        help="azure: South African neural voice (needs a key); edge: the same voice through "
        "the unofficial edge-tts package (no key); windows: an installed offline voice, e.g. "
        "--voice 'Microsoft David Desktop' (placeholder, not South African)",
    )
    a = ap.parse_args()
    sc = scenes(SCRIPT)
    words = sum(len(p.split()) for _, ps in sc for p in ps)
    print(f"{len(sc)} scenes, {words} words (~{words / 150:.1f} min at 150 words/min)")
    if a.dry_run:
        for title, ps in sc:
            print(f"\n# {title}\n{ssml(ps, a.voice, a.rate)}")
        return
    if a.engine == "edge":
        OUT.mkdir(parents=True, exist_ok=True)
        for i, (title, ps) in enumerate(sc, start=1):
            mp3 = OUT / f"scene_{i:02d}.mp3"
            synthesize_edge(plain_text(ps), a.voice, a.rate, mp3)
            print(f"  scene {i:02d} ({title}): {mp3.stat().st_size / 1024:.0f} kB ({a.voice})")
        print(f"wrote {OUT.relative_to(ROOT)}")
        return
    if a.engine == "windows":
        ffmpeg = shutil.which("ffmpeg")
        if ffmpeg is None:
            raise SystemExit("ffmpeg not found")
        voice = a.voice if a.voice != "en-ZA-LukeNeural" else "Microsoft David Desktop"
        OUT.mkdir(parents=True, exist_ok=True)
        for i, (title, ps) in enumerate(sc, start=1):
            wav = OUT / f"scene_{i:02d}.wav"
            synthesize_windows(plain_text(ps), voice, wav)
            mp3 = OUT / f"scene_{i:02d}.mp3"
            subprocess.run(
                [
                    ffmpeg,
                    "-v",
                    "error",
                    "-y",
                    "-i",
                    str(wav),
                    "-ar",
                    "24000",
                    "-ac",
                    "1",
                    "-b:a",
                    "96k",
                    str(mp3),
                ],
                check=True,
            )
            wav.unlink()
            print(f"  scene {i:02d} ({title}): {mp3.stat().st_size / 1024:.0f} kB ({voice})")
        print(f"wrote {OUT.relative_to(ROOT)} (offline voice, placeholder)")
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
