# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

import biotite.structure as struc
import numpy as np
import pytest


@pytest.fixture
def atom_list():
    chain_id = ["A","A","B","B","B"]
    res_id = [1,1,1,1,2]
    res_name = ["ALA","ALA","PRO","PRO","MSE"]
    hetero = [False, False, False, False, True]
    atom_name = ["N", "CA", "O", "CA", "SE"]
    element = ["N","C","O","C","SE"]
    atom_list = []
    for i in range(5):
        atom_list.append(struc.Atom([i,i,i],
                         chain_id = chain_id[i],
                         res_id = res_id[i],
                         res_name = res_name[i],
                         hetero = hetero[i],
                         atom_name = atom_name[i],
                         element = element[i]))
    return atom_list

@pytest.fixture
def atom(atom_list):
    return atom_list[2]

@pytest.fixture
def array(atom_list):
    return struc.array(atom_list)

@pytest.fixture
def stack(array):
    return struc.stack([array, array.copy(), array.copy()])

def test_access(array):
    chain_id = ["A","A","B","B","B"]
    assert array.coord.shape == (5,3)
    assert array.chain_id.tolist() == chain_id
    assert array.get_annotation("chain_id").tolist() == chain_id
    array.add_annotation("test1", dtype=int)
    assert array.test1.tolist() == [0,0,0,0,0]
    with pytest.raises(IndexError):
        array.set_annotation("test2", np.array([0,1,2,3]))

def test_modification(atom, array, stack):
    new_atom = atom
    new_atom.chain_id = "C"
    del array[2]
    assert array.chain_id.tolist() == ["A","A","B","B"]
    array[-1] = new_atom
    assert array.chain_id.tolist() == ["A","A","B","C"]
    del stack[1]
    assert stack.stack_depth() == 2

def test_array_indexing(atom, array):
    filtered_array = array[array.chain_id == "B"]
    assert filtered_array.res_name.tolist() == ["PRO","PRO","MSE"] 
    assert atom == filtered_array[0]
    filtered_array = array[[0,2,4]]
    assert filtered_array.element.tolist() == ["N","O","SE"]

def test_stack_indexing(stack):
    with pytest.raises(IndexError):
        stack[5]
    filtered_stack = stack[0]
    assert type(filtered_stack) == struc.AtomArray
    filtered_stack = stack[0:2, stack.res_name == "PRO"]
    assert filtered_stack.atom_name.tolist() == ["O","CA"]
    filtered_stack = stack[np.array([True,False,True])]
    assert filtered_stack.stack_depth() == 2
    assert filtered_stack.array_length() == 5
    filtered_stack = stack[:,0]
    assert filtered_stack.stack_depth() == 3
    assert filtered_stack.array_length() == 1
    

def test_concatenation(array, stack):
    concat_array = array[2:] + array[:2]
    assert concat_array.chain_id.tolist() == ["B","B","B","A","A"]
    assert concat_array.coord.shape == (5,3)
    concat_stack = stack[:,2:] + stack[:,:2]
    assert concat_array.chain_id.tolist() == ["B","B","B","A","A"]
    assert concat_stack.coord.shape == (3,5,3)

def test_comparison(array):
    mod_array = array.copy()
    assert mod_array == array
    mod_array.coord += 1
    assert mod_array != array
    mod_array = array.copy()
    mod_array.res_name[0] = "UNK"
    mod_array != array

def test_bonds(array):
    assert array.bonds is None
    with pytest.raises(TypeError):
        # Expect a BondList
        array.bonds = 42
    with pytest.raises(IndexError):
        # Expect a BondList with array length as atom count
        array.bonds = struc.BondList(13)
    array.bonds = struc.BondList(5, np.array([(0,1),(0,2),(2,3),(2,4)]))
    assert array.bonds.as_array().tolist() == [[0, 1, 0],
                                               [0, 2, 0],
                                               [2, 3, 0],
                                               [2, 4, 0],]
    filtered_array = array[array.chain_id == "B"]
    assert filtered_array.bonds.as_array().tolist() == [[0, 1, 0],
                                                        [0, 2, 0]]
    concat_array = array + array
    assert concat_array.bonds.as_array().tolist() == [[0, 1, 0],
                                                      [0, 2, 0],
                                                      [2, 3, 0],
                                                      [2, 4, 0],
                                                      [5, 6, 0],
                                                      [5, 7, 0],
                                                      [7, 8, 0],
                                                      [7, 9, 0]]