

# regex patterns
CAMEL = r'^([a-z\d]+)([A-Z\d]+[a-z\d]+)*([A-Z\d]?)?$'
PASCAL = r'^(?<!=[a-z\d])([A-Z\d]+[a-z\d]+)*([A-Z\d]?)?$'
SNAKE = r'^([a-z\d]+)?(_[a-z\d]+)*$'
SNAKE_CAPS = r'^([A-Z\d]+)?(_[A-Z\d]+)*$'