"""
Create on 2-17-2016
@authors: heminggs
"""


class DLListNode:
    """
    Node structure of the linked list.
    Each node has a data field and a next field
    """
    def __init__(self, value):
        self.data = value
        self.prev = None
        self.next = None

    def __str__(self):
        return str(self.data)


class DLIterator:
    def __init__(self, dl_list):
        self._cur = dl_list.head

    def __next__(self):
        if self._cur is None:
            raise StopIteration
        else:
            val = self._cur.data
            self._cur = self._cur.next
            return val


class DLList:
    def __init__(self, data=None):
        self.head = None
        self.tail = None
        self.count = 0
        if data is not None:
            for val in data:
                self.append(val)

    def __len__(self):
        return self.count

    def __str__(self):
        string = "["
        head = self.head
        while head is not None:
            if type(head.data) == str:
                string += "'" + head.data + "'"
            else:
                string += str(head.data)
            head = head.next
            if head is not None:
                string += ", "
        return string + "]"

    def __iter__(self):
        """Return an interable object that can move through the linked list"""
        return DLIterator(self)

    def append(self, value):
        new_node = DLListNode(value)
        if self.count == 0:
            self.head = new_node
            self.tail = new_node
            self.count = 1
        else:
            self.tail.next = new_node
            new_node.prev = self.tail
            self.tail = new_node
            self.count += 1

    def pop(self):
        """Remove and return the last item in the list"""
        val = None
        if self.head is not None:
            val = self.tail.data
            self.count -= 1
            if self.count == 1:
                self.head = None
                self.tail = None
            elif self.count > 1:
                self.tail = self.tail.prev
                self.tail.next = None
        return val

    def insert(self, value, node_after):
        pass

    def remove(self, node):
        pass

    def find(self, value):
        pass

    def __getitem__(self, item):
        pass

    def __setitem__(self, key, value):
        pass

    def reverse(self):
        pass


if __name__ == '__main__':
    # create a linked list with data
    my_list = DLList([1, 2, 3, 4])
    print("start")
    print(my_list)
    print(my_list.head.next)
    my_list.append(5)
    my_list.append([5, 6, 7])
    my_list.append((8, 'a', 10))
    print(my_list)
    print(my_list.pop())
    my_list.append("6")
    print(my_list)
    my_list = DLList()
    print(len(my_list))
    for a in my_list:
        for b in my_list:
            print(a, b)
