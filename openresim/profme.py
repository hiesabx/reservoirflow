import cProfile
import pstats
import io
import re
from functools import wraps
import pandas as pd
import time


def profile(sort_stats="tottime", print_output=True, save_output=True):
    """cProfile decorator.

    Parameters
    ----------
    sort_stats : str, optional, by default "tottime"
        sort table values by ['tottime', 'cumtime', 'ncalls', 'percall']
    print_output : bool, optional, by default True
        print output as a dataframe.
    save_output : bool, optional, by default False
        save output as csv file in the current working dir.

    Backup:
    - Other ways to store csv:
        1. using csv module: result must be a list
                with open(
                    f"profile_{func.__qualname__}.csv", mode="w", newline=""
                ) as f:
                    writer = csv.writer(f, delimiter=",")
                    writer.writerows([r.split(",") for r in result])
        2. using pandas:
                data = io.StringIO("\n".join(result))
                df = pd.read_csv(data, sep=",")
                df.to_csv(f"profile_{func.__qualname__}.csv", index=False)
    """

    def decorator(func):
        @wraps(func)
        def wrapped_func(*args, **kwargs):
            pr = cProfile.Profile()
            s = io.StringIO()
            pr.enable()
            _func = func(*args, **kwargs)
            pr.disable()
            stats = pstats.Stats(pr, stream=s).sort_stats(sort_stats)
            stats.print_stats()
            result = s.getvalue()
            result = [l.strip() for l in re.split("\n+", result)][2:-1]
            result = [re.sub("\s+", ",", r, count=5) for r in result]
            if print_output:
                data = io.StringIO("\n".join(result))
                df = pd.read_csv(data, sep=",")
                print(df)
            if save_output:
                with open(f"profile_{func.__qualname__}.csv", mode="w") as f:
                    f.write("\n".join(result))
            return _func

        return wrapped_func

    return decorator


def timeit(f):
    # reference: https://stackoverflow.com/a/27737385/11549398
    @wraps(f)
    def wrap(*args, **kw):
        ts = time.time()
        result = f(*args, **kw)
        te = time.time()
        # print("func:%r args:[%r, %r] took: %2.4f sec" % (f.__name__, args, kw, te - ts))
        print("func:%r took: %2.4f sec" % (f.__name__, te - ts))
        return result

    return wrap
