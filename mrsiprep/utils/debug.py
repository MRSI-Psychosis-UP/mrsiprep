from __future__ import annotations

import threading
import time
from rich.console import Console
from rich.markup import escape
from rich.theme import Theme


class Debug:
    def __init__(self, verbose: bool = True):
        custom_theme = Theme(
            {
                "success": "green",
                "error": "red",
                "warning": "yellow",
                "failure": "bold red",
                "info": "blue",
                "proc": "violet",
                "debug": "magenta",
            }
        )
        self.console = Console(theme=custom_theme)
        self.verbose = verbose

    def _prepare_message(self, *messages):
        if not messages:
            return "", ""
        prefix = ""
        items = list(messages)
        first = str(items[0])
        if first.strip() == "":
            prefix = first
            items = items[1:]
        text = " ".join(str(message) for message in items) if items else ""
        return prefix, text

    def success(self, *messages):
        if self.verbose:
            prefix, message = self._prepare_message(*messages)
            self.console.print(f"{prefix}[success][ SUCCESS ][/success] {escape(message)}")

    def error(self, *messages):
        if self.verbose:
            prefix, message = self._prepare_message(*messages)
            self.console.print(f"{prefix}[error][  ERROR  ][/error] {escape(message)}")

    def warning(self, *messages):
        if self.verbose:
            prefix, message = self._prepare_message(*messages)
            self.console.print(f"{prefix}[warning][ WARNING ][/warning] {escape(message)}")

    def failure(self, *messages):
        if self.verbose:
            prefix, message = self._prepare_message(*messages)
            self.console.print(f"{prefix}[failure][ FAILURE ][/failure] {escape(message)}")

    def info(self, *messages):
        if self.verbose:
            prefix, message = self._prepare_message(*messages)
            self.console.print(f"{prefix}[info][   INFO  ][/info] {escape(message)}")

    def proc(self, *messages):
        if self.verbose:
            prefix, message = self._prepare_message(*messages)
            self.console.print()
            self.console.print(f"{prefix}[proc][  PROC  ][/proc] {escape(message)}")

    def separator(self):
        if self.verbose:
            self.console.rule()

    def title(self, title: str):
        if self.verbose:
            self.console.rule(title, style="debug")

    def __progress_bar_task(self, duration: int):
        from rich.progress import Progress

        with Progress() as progress:
            task = progress.add_task("[green]Registration...", total=duration)
            for _ in range(duration):
                time.sleep(1)
                progress.update(task, advance=1)

    def run_progress_bar_in_background(self, duration: int):
        thread = threading.Thread(target=self.__progress_bar_task, args=(duration,), daemon=True)
        thread.start()
