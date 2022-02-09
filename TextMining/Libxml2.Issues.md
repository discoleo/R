
# Issues with Libxml2

## Count()

* Valid XPATH:
  xpath1/node1b[count(xpath2) = val]/xpath3;
* Non-functional:
  xpath1/count(xpath2);

The 2nd example is non-functional, although the base functionality already exists!

### Solution
Digging through the libxml2 code:

* Enum xmlXPathObjectType:
  - may need a new type beyond XPATH_NUMBER;
```
Enum xmlXPathObjectType {
    XPATH_UNDEFINED = 0
    XPATH_NODESET = 1
    XPATH_BOOLEAN = 2
    XPATH_NUMBER = 3
    XPATH_STRING = 4
    XPATH_POINT = 5
    XPATH_RANGE = 6
    XPATH_LOCATIONSET = 7
    XPATH_USERS = 8
    XPATH_XSLT_TREE = 9 : An XSLT value tree, non modifiable
}
```

* Code: xmlXPathCacheNewFloat() in xpath.c
  https://github.com/GNOME/libxml2/blob/master/xpath.c

