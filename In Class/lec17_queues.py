from lec15_dllist import DLList

__author__ = 'Kiesligs'

class Queue:
    def __init__(self):
        self.items = DLList

    def __str__(self, items):
        return str(self.items)

    def is_empty(self,items):
        return len(self.items) == 0

    def enqueue(self, item):
        self.items.append(item)

    def dequeue(self):
        if len(self.items) > 0:
            data = self.items.head.data
            self.items.remove(self.item.head)
            return data
        else:
            return None

    def front(self):
        if len(self.items) > 0:
            return self.items.head.data
        else:
            return None

    def size(self):
        return len(self.items)

if __name__ == '__main__':
    q = Queue()
    q.enqueue(1)
    q.enqueue(2)
    q.enqueue(3)
    print(q.dequeue())
    q.enqueue(4)
    print(q.dequeue())
