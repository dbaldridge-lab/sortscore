"""
Top-level CLI dispatcher for sortscore subcommands.
"""
import sys


def _print_usage() -> None:
    print("Usage: sortscore <score|norm> [args]")
    print("")
    print("Commands:")
    print("  score                Run single-experiment scoring workflow")
    print("  norm                 Run batch normalization workflow")
    print("")
    print("Use 'sortscore <command> --help' for command-specific arguments.")


def _dispatch(command_main, argv_tail):
    original_argv = sys.argv[:]
    try:
        sys.argv = [original_argv[0]] + argv_tail
        command_main()
    finally:
        sys.argv = original_argv


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)

    if not args or args[0] in {"-h", "--help"}:
        _print_usage()
        return

    command = args[0]
    command_args = args[1:]

    if command == "score":
        from sortscore.run_analysis import main as score_main
        _dispatch(score_main, command_args)
        return

    if command == "norm":
        from sortscore.run_batch_analysis import main as normalize_main
        _dispatch(normalize_main, command_args)
        return

    print(f"Unknown command: {command}")
    _print_usage()
    raise SystemExit(2)
