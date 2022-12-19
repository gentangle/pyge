"""Decorator to profile time of computation for the core g(aussian)ent(anglement)"""
import io
from pathlib import Path
import pstats
import cProfile

parent_dir = Path(__file__).parent.resolve()

def profileit(func):
    """
    Profile a function with this decorator

    https://stackoverflow.com/questions/5375624/a-decorator-that-profiles-a-method-call-and-logs-the-profiling-result
    """
    def wrapper(*args, **kwargs):
        datafn = func.__name__ + ".profile" # Name the data file sensibly
        prof = cProfile.Profile()
        retval = prof.runcall(func, *args, **kwargs)
        s = io.StringIO()
        ps = pstats.Stats(prof, stream=s).sort_stats('tottime')
        ps.print_stats()
        with open(parent_dir / "profiles" / datafn, 'w+', encoding='utf-8') as out:
            out.write(s.getvalue())
        return retval

    return wrapper
