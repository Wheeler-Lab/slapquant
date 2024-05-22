from queue import Queue
from typing import Generator, Any, Callable
from threading import Thread

def queue_iterator[T](queue: Queue[T]) -> Generator[T, Any, None]:
    counter = 0
    while True:
        counter += 1
        value: T | None = queue.get()
        if value is None:
            break
        yield value

class QueueTee[T]:
    def __init__(self, queue: Queue[T], copies=2, maxsize=None):
        self.input = queue
        self.outputs: list[Queue[T]] = [Queue[T](maxsize=queue.maxsize if maxsize is None else maxsize) for _ in range(copies)]
    
        self.distributor_thread = Thread(target=self.distributor)
        self.distributor_thread.daemon = True
        self.distributor_thread.start()
    
    def distributor(self):
        counter = 0
        while True:
            counter += 1
            item = self.input.get()
            if item is None:
                break
            for output in self.outputs:
                output.put(item)
        for output in self.outputs:
            output.put(None)

class QueueConsumer[T,O](Thread):
    def __init__(self, target: Callable[[T],O], input: Queue[T], output: Queue[O] | None = None, progress=None):
        Thread.__init__(self)
        self.input = input
        self.output = output
        self.target = target
        self.daemon = True

        if progress is not None:
            self._progress = progress
        else:
            self._progress = lambda x: x
    
    def run(self):
        if self.output is not None:
            for item in self._progress(queue_iterator(self.input)):
                processed = self.target(item)
                if processed is not None:
                    self.output.put(processed)
            self.output.put(None)
        else:
            for item in self._progress(queue_iterator(self.input)):
                self.target(item)