"""A Python 3  program to sort an array
according to the order defined by
another array"""

"""A Binary Search based function to find
index of FIRST occurrence of x in arr[].
If x is not present, then it returns -1 """

def first(arr, low, high, x, n) :
    if (high >= low) :
        mid = low + (high - low) // 2;  # (low + high)/2;
        if ((mid == 0 or x > arr[mid-1]) and arr[mid] == x) :
            return mid
        if (x > arr[mid]) :
            return first(arr, (mid + 1), high, x, n)
        return first(arr, low, (mid -1), x, n)

    return -1

# Sort A1[0..m-1] according to the order
# defined by A2[0..n-1].
def sortAccording(A1, A2):
    """The temp array is used to store a copy
    of A1[] and visited[] is used mark the
    visited elements in temp[]."""
    m = len(A1)
    n = len(A2)

    temp = [0] * m
    visited = [0] * m

    for i in range(0, m) :
        temp[i] = A1[i]
        visited[i] = 0

    # Sort elements in temp
    temp.sort()

    # for index of output which is sorted A1[]
    ind = 0

    """Consider all elements of A2[], find
    them in temp[] and copy to A1[] in order."""
    for i in range(0,n) :

        # Find index of the first occurrence
        # of A2[i] in temp
        f = first(temp, 0, m-1, A2[i], m)

        # If not present, no need to proceed
        if (f == -1) :
            continue

        # Copy all occurrences of A2[i] to A1[]
        j = f
        while (j<m and temp[j]==A2[i]) :
            A1[ind] = temp[j];
            ind=ind+1
            visited[j] = 1
            j = j + 1

    # Now copy all items of temp[] which are
    # not present in A2[]
    for i in range(0, m) :
        if (visited[i] == 0) :
            A1[ind] = temp[i]
            ind = ind + 1

# Utility function to print an array
def printArray(arr, n) :
    for i in range(0, n) :
        print(arr[i], end = " ")
    print("")


## Driver program to test above function.
#A1 = [2, 1, 2, 5, 7, 1, 9, 3, 6, 8, 8]
#A2 = [2, 1, 8, 3]
#m = len(A1)
#n = len(A2)
#print("Sorted array is ")
#sortAccording(A1, A2)
#printArray(A1, m)


# This codee is contributed by Nikita Tiwari.