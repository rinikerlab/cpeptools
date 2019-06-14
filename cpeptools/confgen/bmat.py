import numpy as np

def modify_bound_matrix_3(bound_matrix, new_matrix = None, value = None, func = None, modify = "upper", indices = None, check = True, verbose = False):
    """exchange boundmatrix term from new_matrix only if new_matrix term
    is smaller

    Paramters
    ------------
    new_matrix : np.array
    value : float
        only one of value, new_matrx and func should not be None
    func : function
        returns a value as a function of two params i,j which are the corresponding lower bound and upper bound element
    indices : list of int
        if only subset of the bound_matrix that needs to be modified
    modify : str
        whether is to make upper bound smaller or lower bound bigger
    """
    if indices :
        reorder = np.argsort(indices)
        indices.sort()#FIXME  Assumes indices for the ring atom members are in ascending order
        if new_matrix is not None:
            new_matrix = new_matrix[ [[i] for i in reorder] , [reorder]]
        output = np.copy(bound_matrix) #a copy of the full matrix
        bound_matrix = bound_matrix[ [[i] for i in indices], [indices] ]
    counter = 0
    magnitude_of_change = 0
    if new_matrix is None and value is None and func is None:
        return "No valid modification can be done on bound matrix"
    if new_matrix is not None and bound_matrix.shape != new_matrix.shape:
        return "Incompatible shape"
    if not modify in ["upper", "lower"]:
        return "invalid param for modify keyword"

    if modify is "upper" : # lower half is modified only if upper half becomes smaller than it
        for i in range(bound_matrix.shape[0] - 1):
            for j in range(i + 1, bound_matrix.shape[1]):
                if new_matrix is not None:

                    new_matrix[i,j] -= (bound_matrix[i,j]/bound_matrix[j,i])*0.02
                    if check and  bound_matrix[i,j] > new_matrix[i, j]: #Upper bound part



                        counter += 1
                        magnitude_of_change += abs(bound_matrix[i,j] - new_matrix[i, j])
                        bound_matrix[i,j] = new_matrix[i, j]

                        if bound_matrix[j,i] > new_matrix[j, i]: #Lower bound part
                            magnitude_of_change += abs(bound_matrix[j, i] - new_matrix[j, i])
                            counter += 1
                            bound_matrix[j,i] = new_matrix[j, i]
                    elif not check:
                        counter += 1
                        magnitude_of_change += abs(bound_matrix[i,j] - new_matrix[i, j])
                        bound_matrix[i,j] = new_matrix[i, j]


                elif value is not None:
                    if check and bound_matrix[i,j] > value:
                        magnitude_of_change += abs(bound_matrix[i, j] - value)
                        counter += 1
                        bound_matrix[i,j] = value

                        if bound_matrix[j, i] > value:
                            magnitude_of_change += abs(bound_matrix[j, i] - value)
                            counter += 1
                            bound_matrix[j,i] = value

                    elif not check:
                        magnitude_of_change += abs(bound_matrix[i, j] - value)
                        counter += 1
                        bound_matrix[i,j] = value

                elif func is not None:
                    val_min, val_max = bound_matrix[j,i], bound_matrix[i, j]
                    value = func(val_min, val_max)
                    if check and bound_matrix[i,j] > value:
                        magnitude_of_change += abs(bound_matrix[i, j] - value)
                        counter += 1
                        bound_matrix[i,j] = value

                        if bound_matrix[j, i] > value:
                            magnitude_of_change += abs(bound_matrix[j, i] - value)
                            counter += 1
                            bound_matrix[j,i] = value

                    elif not check:
                        magnitude_of_change += abs(bound_matrix[i, j] - value)
                        counter += 1
                        bound_matrix[i,j] = value
                # print(bound_matrix)
    elif modify is "lower" : # upper half is modified only if lower half becomes bigger than it
        for i in range(bound_matrix.shape[0] - 1):
            for j in range(i + 1, bound_matrix.shape[1]):
        # for i in range(1 , bound_matrix.shape[0]):
            # for j in range(i):
                if  new_matrix is not None:
                    if check and bound_matrix[j,i] < new_matrix[j, i]: #Lower bound part
                        magnitude_of_change += abs(bound_matrix[j, i] - new_matrix[j, i])
                        counter += 1
                        bound_matrix[j,i] = new_matrix[j, i]

                        if bound_matrix[i,j] < new_matrix[i, j]: #Upper bound part
                            counter += 1
                            magnitude_of_change += abs(bound_matrix[i,j] - new_matrix[i, j])
                            bound_matrix[i,j] = new_matrix[i, j]
                    elif not check:
                        magnitude_of_change += abs(bound_matrix[j, i] - new_matrix[j, i])
                        counter += 1
                        bound_matrix[j,i] = new_matrix[j, i]

                elif value is not None:
                    if check and bound_matrix[j,i] < value:
                        magnitude_of_change += abs(bound_matrix[j, i] - value)
                        counter += 1
                        bound_matrix[j,i] = value

                        if bound_matrix[i, j] < value:
                            magnitude_of_change += abs(bound_matrix[i, j] - value)
                            counter += 1
                            bound_matrix[i,j] = value

                    elif not check:
                        magnitude_of_change += abs(bound_matrix[j, i] - value)
                        counter += 1
                        bound_matrix[j,i] = value

                elif func is not None:
                    val_min, val_max = bound_matrix[j,i], bound_matrix[i, j]
                    value = func(val_min, val_max)
                    if check and bound_matrix[j,i] < value:
                        magnitude_of_change += abs(bound_matrix[j, i] - value)
                        counter += 1
                        bound_matrix[j,i] = value

                        # value = func(i, j)
                        if bound_matrix[i, j] < value:
                            magnitude_of_change += abs(bound_matrix[i, j] - value)
                            counter += 1
                            bound_matrix[i,j] = value

                    elif not check:
                        magnitude_of_change += abs(bound_matrix[j, i] - value)
                        counter += 1
                        bound_matrix[j,i] = value

    if verbose:
        print("Number of changes : {}".format(counter))
        print("Amount of changes : {}".format(magnitude_of_change))
    if indices :
        if verbose :
            print(indices, bound_matrix)
        output[ [[i] for i in indices], [indices] ] = bound_matrix
        bound_matrix = output
    return bound_matrix

# def modify_bound_matrix(bound_matrix, new_matrix = None, value = None, func = None, modify = "upper", indices = None, check = True, verbose = False):
#     """exchange boundmatrix term from new_matrix only if new_matrix term
#     is smaller
#
#     Paramters
#     ------------
#     new_matrix : np.array
#     value : float
#         only one of value, new_matrx and func should not be None
#     func : function
#         returns a value as a function of two params i,j which are the corresponding lower bound and upper bound element
#     indices : list of int
#         if only subset of the bound_matrix that needs to be modified
#     modify : str
#         whether is to make upper bound smaller or lower bound bigger
#     """
#     if indices :
#         reorder = np.argsort(indices)
#         indices.sort()#FIXME  Assumes indices for the ring atom members are in ascending order
#         if new_matrix is not None:
#             new_matrix = new_matrix[ [[i] for i in reorder] , [reorder]]
#         output = np.copy(bound_matrix) #a copy of the full matrix
#         bound_matrix = bound_matrix[ [[i] for i in indices], [indices] ]
#     counter = 0
#     magnitude_of_change = 0
#     if new_matrix is None and value is None and func is None:
#         return "No valid modification can be done on bound matrix"
#     if new_matrix is not None and bound_matrix.shape != new_matrix.shape:
#         return "Incompatible shape"
#     if not modify in ["upper", "lower"]:
#         return "invalid param for modify keyword"
#
#     if modify is "upper" : # lower half is modified only if upper half becomes smaller than it
#         for i in range(bound_matrix.shape[0] - 1):
#             for j in range(i + 1, bound_matrix.shape[1]):
#                 if new_matrix is not None:
#                     if check and  bound_matrix[i,j] > new_matrix[i, j]: #Upper bound part
#                         counter += 1
#                         magnitude_of_change += abs(bound_matrix[i,j] - new_matrix[i, j])
#                         bound_matrix[i,j] = new_matrix[i, j]
#
#                         if bound_matrix[j,i] > new_matrix[j, i]: #Lower bound part
#                             magnitude_of_change += abs(bound_matrix[j, i] - new_matrix[j, i])
#                             counter += 1
#                             bound_matrix[j,i] = new_matrix[j, i]
#                     elif not check:
#                         counter += 1
#                         magnitude_of_change += abs(bound_matrix[i,j] - new_matrix[i, j])
#                         bound_matrix[i,j] = new_matrix[i, j]
#
#
#                 elif value is not None:
#                     if check and bound_matrix[i,j] > value:
#                         magnitude_of_change += abs(bound_matrix[i, j] - value)
#                         counter += 1
#                         bound_matrix[i,j] = value
#
#                         if bound_matrix[j, i] > value:
#                             magnitude_of_change += abs(bound_matrix[j, i] - value)
#                             counter += 1
#                             bound_matrix[j,i] = value
#
#                     elif not check:
#                         magnitude_of_change += abs(bound_matrix[i, j] - value)
#                         counter += 1
#                         bound_matrix[i,j] = value
#
#                 elif func is not None:
#                     val_min, val_max = bound_matrix[j,i], bound_matrix[i, j]
#                     value = func(val_min, val_max)
#                     if check and bound_matrix[i,j] > value:
#                         magnitude_of_change += abs(bound_matrix[i, j] - value)
#                         counter += 1
#                         bound_matrix[i,j] = value
#
#                         if bound_matrix[j, i] > value:
#                             magnitude_of_change += abs(bound_matrix[j, i] - value)
#                             counter += 1
#                             bound_matrix[j,i] = value
#
#                     elif not check:
#                         magnitude_of_change += abs(bound_matrix[i, j] - value)
#                         counter += 1
#                         bound_matrix[i,j] = value
#                 # print(bound_matrix)
#     elif modify is "lower" : # upper half is modified only if lower half becomes bigger than it
#         for i in range(bound_matrix.shape[0] - 1):
#             for j in range(i + 1, bound_matrix.shape[1]):
#         # for i in range(1 , bound_matrix.shape[0]):
#             # for j in range(i):
#                 if  new_matrix is not None:
#                     if check and bound_matrix[j,i] < new_matrix[j, i]: #Lower bound part
#                         magnitude_of_change += abs(bound_matrix[j, i] - new_matrix[j, i])
#                         counter += 1
#                         bound_matrix[j,i] = new_matrix[j, i]
#
#                         if bound_matrix[i,j] < new_matrix[i, j]: #Upper bound part
#                             counter += 1
#                             magnitude_of_change += abs(bound_matrix[i,j] - new_matrix[i, j])
#                             bound_matrix[i,j] = new_matrix[i, j]
#                     elif not check:
#                         magnitude_of_change += abs(bound_matrix[j, i] - new_matrix[j, i])
#                         counter += 1
#                         bound_matrix[j,i] = new_matrix[j, i]
#
#                 elif value is not None:
#                     if check and bound_matrix[j,i] < value:
#                         magnitude_of_change += abs(bound_matrix[j, i] - value)
#                         counter += 1
#                         bound_matrix[j,i] = value
#
#                         if bound_matrix[i, j] < value:
#                             magnitude_of_change += abs(bound_matrix[i, j] - value)
#                             counter += 1
#                             bound_matrix[i,j] = value
#
#                     elif not check:
#                         magnitude_of_change += abs(bound_matrix[j, i] - value)
#                         counter += 1
#                         bound_matrix[j,i] = value
#
#                 elif func is not None:
#                     val_min, val_max = bound_matrix[j,i], bound_matrix[i, j]
#                     value = func(val_min, val_max)
#                     if check and bound_matrix[j,i] < value:
#                         magnitude_of_change += abs(bound_matrix[j, i] - value)
#                         counter += 1
#                         bound_matrix[j,i] = value
#
#                         # value = func(i, j)
#                         if bound_matrix[i, j] < value:
#                             magnitude_of_change += abs(bound_matrix[i, j] - value)
#                             counter += 1
#                             bound_matrix[i,j] = value
#
#                     elif not check:
#                         magnitude_of_change += abs(bound_matrix[j, i] - value)
#                         counter += 1
#                         bound_matrix[j,i] = value
#
#     if verbose:
#         print("Number of changes : {}".format(counter))
#         print("Amount of changes : {}".format(magnitude_of_change))
#     if indices :
#         if verbose :
#             print(indices, bound_matrix)
#         output[ [[i] for i in indices], [indices] ] = bound_matrix
#         bound_matrix = output
#     return bound_matrix
#
# #!!!!!! this is one I am currently experimenting!!!!!
# def modify_bound_matrix(bound_matrix, new_matrix = None, value = None, func = None, modify = "upper", indices = None, check = True, verbose = False):
#     """exchange boundmatrix term from new_matrix only if new_matrix term
#     is smaller
#
#     Paramters
#     ------------
#     new_matrix : np.array
#     value : float
#         only one of value, new_matrx and func should not be None
#     func : function
#         returns a value as a function of two params i,j which are the corresponding lower bound and upper bound element
#     indices : list of int
#         if only subset of the bound_matrix that needs to be modified
#     modify : str
#         whether is to make upper bound smaller or lower bound bigger
#     """
#     if indices :
#         reorder = np.argsort(indices)
#         indices.sort()#FIXME  Assumes indices for the ring atom members are in ascending order
#         if new_matrix is not None:
#             new_matrix = new_matrix[ [[i] for i in reorder] , [reorder]]
#         output = np.copy(bound_matrix) #a copy of the full matrix
#         bound_matrix = bound_matrix[ [[i] for i in indices], [indices] ]
#     counter = 0
#     magnitude_of_change = 0
#     if new_matrix is None and value is None and func is None:
#         return "No valid modification can be done on bound matrix"
#     if new_matrix is not None and bound_matrix.shape != new_matrix.shape:
#         return "Incompatible shape"
#     if not modify in ["upper", "lower"]:
#         return "invalid param for modify keyword"
#
#     GEN_DIST_TOL = 0.06
#     if modify is "upper" : # lower half is modified only if upper half becomes smaller than it
#         for i in range(bound_matrix.shape[0] - 1):
#             for j in range(i + 1, bound_matrix.shape[1]):
#                 if new_matrix is not None:
#                     if not check:
#                         counter += 1
#                         magnitude_of_change += abs(bound_matrix[i,j] - new_matrix[i, j])
#                         bound_matrix[i,j] = new_matrix[i, j]
#
#                     elif bound_matrix[j, i] < new_matrix[j, i] and new_matrix[i, j] <= bound_matrix[i, j]:
#                         if new_matrix[j,i] - bound_matrix[j,i] < bound_matrix[i, j] - new_matrix[i , j]:
#                             bound_matrix[j, i] = new_matrix[j, i]
#                         bound_matrix[i, j] = new_matrix[i, j]
#
#
#                     # elif new_matrix[j, i] < bound_matrix[j, i] :
#                     #     bound_matrix[j, i] = new_matrix[j, i]
#                     #
#                     # elif bound_matrix[i, j] < new_matrix[i, j]:
#                     #     bound_matrix[i, j] = new_matrix[i, j]
#
#
#                 elif value is not None:
#                     if check and bound_matrix[i,j] > value:
#                         magnitude_of_change += abs(bound_matrix[i, j] - value)
#                         counter += 1
#                         bound_matrix[i,j] = value
#
#                         if bound_matrix[j, i] > value:
#                             magnitude_of_change += abs(bound_matrix[j, i] - value)
#                             counter += 1
#                             bound_matrix[j,i] = value
#
#                     elif not check:
#                         magnitude_of_change += abs(bound_matrix[i, j] - value)
#                         counter += 1
#                         bound_matrix[i,j] = value
#
#                 elif func is not None:
#                     val_min, val_max = bound_matrix[j,i], bound_matrix[i, j]
#                     value = func(val_min, val_max)
#                     if check and bound_matrix[i,j] > value:
#                         magnitude_of_change += abs(bound_matrix[i, j] - value)
#                         counter += 1
#                         bound_matrix[i,j] = value
#
#                         if bound_matrix[j, i] > value:
#                             magnitude_of_change += abs(bound_matrix[j, i] - value)
#                             counter += 1
#                             bound_matrix[j,i] = value
#
#                     elif not check:
#                         magnitude_of_change += abs(bound_matrix[i, j] - value)
#                         counter += 1
#                         bound_matrix[i,j] = value
#                 # print(bound_matrix)
#     elif modify is "lower" : # upper half is modified only if lower half becomes bigger than it
#         for i in range(bound_matrix.shape[0] - 1):
#             for j in range(i + 1, bound_matrix.shape[1]):
#         # for i in range(1 , bound_matrix.shape[0]):
#             # for j in range(i):
#                 if  new_matrix is not None:
#                     if check and bound_matrix[j,i] < new_matrix[j, i]: #Lower bound part
#                         magnitude_of_change += abs(bound_matrix[j, i] - new_matrix[j, i])
#                         counter += 1
#                         bound_matrix[j,i] = new_matrix[j, i]
#
#                         if bound_matrix[i,j] < new_matrix[i, j]: #Upper bound part
#                             counter += 1
#                             magnitude_of_change += abs(bound_matrix[i,j] - new_matrix[i, j])
#                             bound_matrix[i,j] = new_matrix[i, j]
#                     elif not check:
#                         magnitude_of_change += abs(bound_matrix[j, i] - new_matrix[j, i])
#                         counter += 1
#                         bound_matrix[j,i] = new_matrix[j, i]
#
#                 elif value is not None:
#                     if check and bound_matrix[j,i] < value:
#                         magnitude_of_change += abs(bound_matrix[j, i] - value)
#                         counter += 1
#                         bound_matrix[j,i] = value
#
#                         if bound_matrix[i, j] < value:
#                             magnitude_of_change += abs(bound_matrix[i, j] - value)
#                             counter += 1
#                             bound_matrix[i,j] = value
#
#                     elif not check:
#                         magnitude_of_change += abs(bound_matrix[j, i] - value)
#                         counter += 1
#                         bound_matrix[j,i] = value
#
#                 elif func is not None:
#                     val_min, val_max = bound_matrix[j,i], bound_matrix[i, j]
#                     value = func(val_min, val_max)
#                     if check and bound_matrix[j,i] < value:
#                         magnitude_of_change += abs(bound_matrix[j, i] - value)
#                         counter += 1
#                         bound_matrix[j,i] = value
#
#                         # value = func(i, j)
#                         if bound_matrix[i, j] < value:
#                             magnitude_of_change += abs(bound_matrix[i, j] - value)
#                             counter += 1
#                             bound_matrix[i,j] = value
#
#                     elif not check:
#                         magnitude_of_change += abs(bound_matrix[j, i] - value)
#                         counter += 1
#                         bound_matrix[j,i] = value
#
#     if verbose:
#         print("Number of changes : {}".format(counter))
#         print("Amount of changes : {}".format(magnitude_of_change))
#     if indices :
#         if verbose :
#             print(indices, bound_matrix)
#         output[ [[i] for i in indices], [indices] ] = bound_matrix
#         bound_matrix = output
#     return bound_matrix


# #trying new, commented out a copy of the old
# #this gives very sharp distributions with not so high eccentricity
def modify_bound_matrix_2(bound_matrix, new_matrix = None, value = None, func = None, modify = "upper", indices = None, check = True, verbose = False):
    """exchange boundmatrix term from new_matrix only if new_matrix term
    is smaller

    Paramters
    ------------
    new_matrix : np.array
    value : float
        only one of value, new_matrx and func should not be None
    func : function
        returns a value as a function of two params i,j which are the corresponding lower bound and upper bound element
    indices : list of int
        if only subset of the bound_matrix that needs to be modified
    modify : str
        whether is to make upper bound smaller or lower bound bigger
    """
    GEN_DIST_TOL = 0.06
    if indices :
        reorder = np.argsort(indices)
        indices.sort()#FIXME  Assumes indices for the ring atom members are in ascending order
        if new_matrix is not None:
            new_matrix = new_matrix[ [[i] for i in reorder] , [reorder]]
        output = np.copy(bound_matrix) #a copy of the full matrix
        bound_matrix = bound_matrix[ [[i] for i in indices], [indices] ]
    counter = 0
    magnitude_of_change = 0
    if new_matrix is None and value is None and func is None:
        return "No valid modification can be done on bound matrix"
    if new_matrix is not None and bound_matrix.shape != new_matrix.shape:
        return "Incompatible shape"
    if not modify in ["upper", "lower"]:
        return "invalid param for modify keyword"

    if modify is "upper" : # lower half is modified only if upper half becomes smaller than it
        for i in range(bound_matrix.shape[0] - 1):
            for j in range(i + 1, bound_matrix.shape[1]):
                if new_matrix is not None:
                    if not check:
                        counter += 1
                        magnitude_of_change += abs(bound_matrix[i,j] - new_matrix[i, j])
                        bound_matrix[i,j] = new_matrix[i, j]

                    elif bound_matrix[j, i] < new_matrix[j, i] and new_matrix[i, j] <= bound_matrix[i, j]:
                        bound_matrix[i, j] = new_matrix[i, j]

                    elif new_matrix[j, i] < bound_matrix[j, i] :
                        bound_matrix[i ,j] = bound_matrix[j ,i]
                        bound_matrix[j, i] = new_matrix[j, i]

                    elif bound_matrix[i, j] < new_matrix[i, j]:
                        bound_matrix[j, i] = bound_matrix[i ,j]
                        bound_matrix[i, j] = new_matrix[i, j]

                    elif bound_matrix[j, i] == new_matrix[j, i]: #change upper bound
                        counter += 1
                        magnitude_of_change += abs(bound_matrix[i,j] - new_matrix[i, j])
                        bound_matrix[i,j] = new_matrix[i, j] + GEN_DIST_TOL

                elif value is not None:
                    if check and bound_matrix[i,j] > value:
                        magnitude_of_change += abs(bound_matrix[i, j] - value)
                        counter += 1
                        bound_matrix[i,j] = value

                        if bound_matrix[j, i] > value:
                            magnitude_of_change += abs(bound_matrix[j, i] - value)
                            counter += 1
                            bound_matrix[j,i] = value

                    elif not check:
                        magnitude_of_change += abs(bound_matrix[i, j] - value)
                        counter += 1
                        bound_matrix[i,j] = value

                elif func is not None:
                    val_min, val_max = bound_matrix[j,i], bound_matrix[i, j]
                    value = func(val_min, val_max)
                    if check and bound_matrix[i,j] > value:
                        magnitude_of_change += abs(bound_matrix[i, j] - value)
                        counter += 1
                        bound_matrix[i,j] = value

                        if bound_matrix[j, i] > value:
                            magnitude_of_change += abs(bound_matrix[j, i] - value)
                            counter += 1
                            bound_matrix[j,i] = value

                    elif not check:
                        magnitude_of_change += abs(bound_matrix[i, j] - value)
                        counter += 1
                        bound_matrix[i,j] = value
                # print(bound_matrix)
    elif modify is "lower" : # upper half is modified only if lower half becomes bigger than it
        for i in range(bound_matrix.shape[0] - 1):
            for j in range(i + 1, bound_matrix.shape[1]):
        # for i in range(1 , bound_matrix.shape[0]):
            # for j in range(i):
                if  new_matrix is not None:
                    if check and bound_matrix[j,i] < new_matrix[j, i]: #Lower bound part
                        magnitude_of_change += abs(bound_matrix[j, i] - new_matrix[j, i])
                        counter += 1
                        bound_matrix[j,i] = new_matrix[j, i]

                        if bound_matrix[i,j] < new_matrix[i, j]: #Upper bound part
                            counter += 1
                            magnitude_of_change += abs(bound_matrix[i,j] - new_matrix[i, j])
                            bound_matrix[i,j] = new_matrix[i, j]
                    elif not check:
                        magnitude_of_change += abs(bound_matrix[j, i] - new_matrix[j, i])
                        counter += 1
                        bound_matrix[j,i] = new_matrix[j, i]

                elif value is not None:
                    if check and bound_matrix[j,i] < value:
                        magnitude_of_change += abs(bound_matrix[j, i] - value)
                        counter += 1
                        bound_matrix[j,i] = value

                        if bound_matrix[i, j] < value:
                            magnitude_of_change += abs(bound_matrix[i, j] - value)
                            counter += 1
                            bound_matrix[i,j] = value

                    elif not check:
                        magnitude_of_change += abs(bound_matrix[j, i] - value)
                        counter += 1
                        bound_matrix[j,i] = value

                elif func is not None:
                    val_min, val_max = bound_matrix[j,i], bound_matrix[i, j]
                    value = func(val_min, val_max)
                    if check and bound_matrix[j,i] < value:
                        magnitude_of_change += abs(bound_matrix[j, i] - value)
                        counter += 1
                        bound_matrix[j,i] = value

                        # value = func(i, j)
                        if bound_matrix[i, j] < value:
                            magnitude_of_change += abs(bound_matrix[i, j] - value)
                            counter += 1
                            bound_matrix[i,j] = value

                    elif not check:
                        magnitude_of_change += abs(bound_matrix[j, i] - value)
                        counter += 1
                        bound_matrix[j,i] = value

    if verbose:
        print("Number of changes : {}".format(counter))
        print("Amount of changes : {}".format(magnitude_of_change))
    if indices :
        if verbose :
            print(indices, bound_matrix)
        output[ [[i] for i in indices], [indices] ] = bound_matrix
        bound_matrix = output
    return bound_matrix


def modify_bound_matrix(bound_matrix, new_matrix = None, value = None, func = None, modify = "upper", indices = None, check = True, verbose = False):
    """exchange boundmatrix term from new_matrix only if new_matrix term
    is smaller

    Paramters
    ------------
    new_matrix : np.array
    value : float
        only one of value, new_matrx and func should not be None
    func : function
        returns a value as a function of two params i,j which are the corresponding lower bound and upper bound element
    indices : list of int
        if only subset of the bound_matrix that needs to be modified
    modify : str
        whether is to make upper bound smaller or lower bound bigger
    """
    if indices :
        reorder = np.argsort(indices)
        indices.sort()#FIXME  Assumes indices for the ring atom members are in ascending order
        if new_matrix is not None:
            new_matrix = new_matrix[ [[i] for i in reorder] , [reorder]]
        output = np.copy(bound_matrix) #a copy of the full matrix
        bound_matrix = bound_matrix[ [[i] for i in indices], [indices] ]
    counter = 0
    magnitude_of_change = 0
    if new_matrix is None and value is None and func is None:
        return "No valid modification can be done on bound matrix"
    if new_matrix is not None and bound_matrix.shape != new_matrix.shape:
        return "Incompatible shape"
    if not modify in ["upper", "lower"]:
        return "invalid param for modify keyword"

    if modify is "upper" : # lower half is modified only if upper half becomes smaller than it
        for i in range(bound_matrix.shape[0] - 1):
            for j in range(i + 1, bound_matrix.shape[1]):
                if new_matrix is not None:
                    if check and  bound_matrix[i,j] > new_matrix[i, j]: #Upper bound part
                        counter += 1
                        magnitude_of_change += abs(bound_matrix[i,j] - new_matrix[i, j])
                        bound_matrix[i,j] = new_matrix[i, j]

                        if bound_matrix[j,i] > new_matrix[j, i]: #Lower bound part
                            magnitude_of_change += abs(bound_matrix[j, i] - new_matrix[j, i])
                            counter += 1
                            bound_matrix[j,i] = new_matrix[j, i]
                    elif not check:
                        counter += 1
                        magnitude_of_change += abs(bound_matrix[i,j] - new_matrix[i, j])
                        bound_matrix[i,j] = new_matrix[i, j]


                elif value is not None:
                    if check and bound_matrix[i,j] > value:
                        magnitude_of_change += abs(bound_matrix[i, j] - value)
                        counter += 1
                        bound_matrix[i,j] = value

                        if bound_matrix[j, i] > value:
                            magnitude_of_change += abs(bound_matrix[j, i] - value)
                            counter += 1
                            bound_matrix[j,i] = value

                    elif not check:
                        magnitude_of_change += abs(bound_matrix[i, j] - value)
                        counter += 1
                        bound_matrix[i,j] = value

                elif func is not None:
                    val_min, val_max = bound_matrix[j,i], bound_matrix[i, j]
                    value = func(val_min, val_max)
                    if check and bound_matrix[i,j] > value:
                        magnitude_of_change += abs(bound_matrix[i, j] - value)
                        counter += 1
                        bound_matrix[i,j] = value

                        if bound_matrix[j, i] > value:
                            magnitude_of_change += abs(bound_matrix[j, i] - value)
                            counter += 1
                            bound_matrix[j,i] = value

                    elif not check:
                        magnitude_of_change += abs(bound_matrix[i, j] - value)
                        counter += 1
                        bound_matrix[i,j] = value
                # print(bound_matrix)
    elif modify is "lower" : # upper half is modified only if lower half becomes bigger than it
        for i in range(bound_matrix.shape[0] - 1):
            for j in range(i + 1, bound_matrix.shape[1]):
        # for i in range(1 , bound_matrix.shape[0]):
            # for j in range(i):
                if  new_matrix is not None:
                    if check and bound_matrix[j,i] < new_matrix[j, i]: #Lower bound part
                        magnitude_of_change += abs(bound_matrix[j, i] - new_matrix[j, i])
                        counter += 1
                        bound_matrix[j,i] = new_matrix[j, i]

                        if bound_matrix[i,j] < new_matrix[i, j]: #Upper bound part
                            counter += 1
                            magnitude_of_change += abs(bound_matrix[i,j] - new_matrix[i, j])
                            bound_matrix[i,j] = new_matrix[i, j]
                    elif not check:
                        magnitude_of_change += abs(bound_matrix[j, i] - new_matrix[j, i])
                        counter += 1
                        bound_matrix[j,i] = new_matrix[j, i]

                elif value is not None:
                    if check and bound_matrix[j,i] < value:
                        magnitude_of_change += abs(bound_matrix[j, i] - value)
                        counter += 1
                        bound_matrix[j,i] = value

                        if bound_matrix[i, j] < value:
                            magnitude_of_change += abs(bound_matrix[i, j] - value)
                            counter += 1
                            bound_matrix[i,j] = value

                    elif not check:
                        magnitude_of_change += abs(bound_matrix[j, i] - value)
                        counter += 1
                        bound_matrix[j,i] = value

                elif func is not None:
                    val_min, val_max = bound_matrix[j,i], bound_matrix[i, j]
                    value = func(val_min, val_max)
                    if check and bound_matrix[j,i] < value:
                        magnitude_of_change += abs(bound_matrix[j, i] - value)
                        counter += 1
                        bound_matrix[j,i] = value

                        # value = func(i, j)
                        if bound_matrix[i, j] < value:
                            magnitude_of_change += abs(bound_matrix[i, j] - value)
                            counter += 1
                            bound_matrix[i,j] = value

                    elif not check:
                        magnitude_of_change += abs(bound_matrix[j, i] - value)
                        counter += 1
                        bound_matrix[j,i] = value

    if verbose:
        print("Number of changes : {}".format(counter))
        print("Amount of changes : {}".format(magnitude_of_change))
    if indices :
        if verbose :
            print(indices, bound_matrix)
        output[ [[i] for i in indices], [indices] ] = bound_matrix
        bound_matrix = output
    return bound_matrix

def check_bmat(bmat):
    for i in range(bmat.shape[0] - 1):
        for j in range(i + 1, bmat.shape[1]):
            if bmat[i,j] < bmat[j,i]:
                raise ValueError("Lower bound {} bigger than upper bound {}".format(i,j))
    print("Boundmatrix is fine")

def bmat_from_conformer(conf):
    """
    Parameters
    -----------
    conf RDKit conformer object
    """
    positions = conf.GetPositions()
    num_atoms = len(positions)
    bmat = cdist(positions, positions)
    return bmat


def get_ring_bond_length_list(bmat, ring_indices, scale_factor = 1.0):
    indices = ring_indices + [ring_indices[0]]
    out = []
    for i in range(len(indices) - 1):
        col_idx = min(indices[i], indices[i+1])
        row_idx = max(indices[i], indices[i+1])
        out.append(bmat[row_idx, col_idx] * scale_factor) #min dis
        #out.append(bmat[col_idx, row_idx] * sacle_factor) #max dis
    return out
