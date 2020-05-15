import random
def generate_one_random_p_m(l):
    alphabet = [0,1,2,3]
    p_m = []
    for i in range(l):
        p_m.append(random.choice(alphabet))
    return p_m

def generate_messages(q,l):
    p_ms = []
    for i in range(q):
        p_ms.append(generate_one_random_p_m(l))
    return p_ms
