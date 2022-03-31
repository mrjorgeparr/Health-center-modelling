from fase3 import Map
from fase3 import HealthCenter


def test():
    m = Map()
    centers = ["A", "B", "C", "D", "E", "F", "G"]
    for center in centers:
        m.addHealthCenter(HealthCenter(center))

    m.addConnection(m.centers[0], m.centers[1], 10)  #A,B 10
    m.addConnection(m.centers[0], m.centers[2], 6)  #A,C 6
    m.addConnection(m.centers[0], m.centers[3], 1)  #A,D 1
    m.addConnection(m.centers[1], m.centers[5], 20)  #B,F 20
    m.addConnection(m.centers[2], m.centers[3], 3)  #C,D 3
    m.addConnection(m.centers[2], m.centers[6], 5)  #C,G 5
    m.addConnection(m.centers[2], m.centers[4], 2)  #C,E
    m.addConnection(m.centers[2], m.centers[5], 4)  #C,F 4
    m.addConnection(m.centers[3], m.centers[6], 7)  #D,G 7
    m.addConnection(m.centers[4], m.centers[6], 2)  #E,G 2
    m.addConnection(m.centers[5], m.centers[4], 1)  #F,E 1
    aux = HealthCenter("J")
    print("Bellman-ford test's\n")
    print("Test 1: from a node to itself")
    minimum_path, d = m.minimumPathBF(m.centers[0], m.centers[0])
    #assert len(minimum_path) == 1, "test 1 BF: incorrect path"
    #assert minimum_path[0] == 0, "test 1 BF: incorrect path"
    print("\nPath: ", minimum_path)
    print("\nDistance: ", d)
    print("\n")
    print("Test 2: from center 'A', to center 'E'")
    minimum_path, d = m.minimumPathBF(m.centers[0], m.centers[4])
    #assert len(correct_path) == len(minimum_path), "test 2 BF: incorrect path"
    print("\nPath: ", minimum_path)
    #assert d == 6
    print("\nDistance: ", d)
    print("\n")
    print("Test 3: from a valid center to a non existing goal center")
    minimum_path, d = m.minimumPathBF(m.centers[0], aux)
    print("\nPath: ", minimum_path)
    #assert d == 0, "test 3 BF: incorrect distance"
    print("\nDistance: ", d)
    print("\n")
    print("Test 4: from  a non existing center to a valid goal center")
    minimum_path, d = m.minimumPathBF(aux, m.centers[0])
    #assert len(minimum_path) == 0, "test 4 BF: incorrect path"
    print("\nPath: ", minimum_path)
    #assert d == 0, "test 4 BF: incorrect distance"
    print("\nDistance: ", d)
    print("\n\n")
    print("Floyd-Warshal test's\n")
    print("Test 1: from a node to itself")
    minimum_path, d = m.minimumPathFW(m.centers[0], m.centers[0])
    #assert len(minimum_path) == 1
    #assert minimum_path[0] == 0, "test 1 FW: incorrect path"
    print("\nPath: ", minimum_path)
    assert d == 0, "test 1 FW: incorrect distance"
    print("\nDistance: ", d)
    print("\n")
    print("Test 2: from center 'A', to center 'E'")
    minimum_path, d = m.minimumPathFW(m.centers[0], m.centers[4])
    #correct_path = [0, 3, 2, 4]
    print("\nPath: ", minimum_path)
    print("\nDistance: ", d)
    print("\n")
    print("Test 3: from a valid center to a non existing goal center")
    minimum_path, d = m.minimumPathFW(m.centers[0], aux)
    assert len(minimum_path) == 0, "test 3 FW: invalid path"
    print("\nPath: ", minimum_path)
    assert d == 0, "test 3 FW: invalid distance"
    print("\nDistance: ", d)
    print("\n")
    print("Test 4: from  a non existing center to a valid goal center")
    minimum_path, d = m.minimumPathFW(aux, m.centers[0])
    assert len(minimum_path) == 0, "test 4 FW: incorrect path"
    print("\nPath: ", minimum_path)
    assert d == 0, "test 4 FW: invalid distance"
    print("\nDistance: ", d)
    print("\n")



test()


def tbackwards_conversion(time):
    hours = time // 60
    minutes = time - hours * 60
    if minutes < 10:
        return str(hours) + ":" + "0" + str(minutes)
    if hours < 10:
        if minutes >= 10:
            return "0" + str(hours) + ":" + str(minutes)
    else:
        return str(hours) + ":" + str(minutes)
def tconversion(time):
    t = time.split(":")
    return int(t[0]) * 60 + int(t[1])


a = "09:00"
b = tconversion(a)
c = tbackwards_conversion(b)
print(c)