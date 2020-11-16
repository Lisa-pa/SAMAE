# These functions assume that you have an obj folder in your current working
# directory, which will be used to store the objects.
# Note that pickle.HIGHEST_PROTOCOL is a binary format, which could not be always
# convenient, but is good for performance. Protocol 0 is a text format.​
# These two functions have been developed by Antonio J. Morales-Artacho
def save_obj(obj, path):
    import pickle
    with open(path + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name, path):
    import pickle
    with open(path + name + '.pkl', 'rb') as f:
        return pickle.load(f)