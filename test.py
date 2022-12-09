



class MyClass:
    files: dict[int, str] = {}

    @classmethod
    def run(cls) -> None:
        for i in range(10):
            cls.update(i, chr(i))

        print(cls.files)

    @classmethod
    def update(cls, i: int, char: str) -> None:
        cls.files[i] = char


MyClass.run()
MyClass.run()




