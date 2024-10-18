from collections.abc import Generator
import collections


class Tee:

    def __init__(
        self,
        gen: Generator,
        yld: callable = lambda x: x,
        **kwargs
    ) -> None:

        self.gen = gen
        self.yld = yld
        self.kwargs = kwargs
        self.cached = collections.defaultdict(list)


    def __iter__(self) -> Generator:

        for item in self.gen:

            for name, callback in self.kwargs.items():

                self.cached[name].append(callback(item))

            yield self.yld(item)

        self.cached = dict(self.cached)
