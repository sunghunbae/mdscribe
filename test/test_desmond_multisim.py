import pytest
from mdscribe.desmond import Multisim


def test_Multisim_variable():
    # parse_string returns ['barostat.tau']
    o1 = Multisim.variable.parse_string("barostat.tau")[0]
    o2 = Multisim.variable.parse_string("barostat  .tau")[0]
    assert o1 == o2
    

def test_Multisim_expression_1():
    i = "task {} simulate { meta = [{a=1 b=3 c=[7 8 9]}] f = {} } simulate {n=2} simulate {n = 3}"
    D = Multisim(string=i).to_dot()
    assert D.simulate[0].meta[0].b == str(3)
    assert D.simulate[0].meta[0].c[1] == str(8)
    assert D.simulate[-1].n == str(3)


def test_Multisim_expression_2():
    i = "simulate { effect_if = [[ a ] b ] }"
    D = Multisim(string=i).to_dot()
    assert D.simulate.effect_if[0] == ['a']
    assert D.simulate.effect_if[1] == 'b'
    

def test_Multisim_expression_3():
    i = """ensemble = {
        brownie = {
            delta_max = 0.1
        }
        class = NVT
        method = Brownie
        thermostat = {
            tau = 1.0
        }
    }"""
    D = Multisim(string=i)
    assert D.dot.ensemble.brownie.delta_max == str(0.1)
    # assert D.dot.ensemble.class == 'NVT' # <--- it raises SyntaxError
    assert D.dot.ensemble.method == 'Brownie'
    assert D.dot.ensemble.thermostat.tau == str(1.0)
    D.write()


def test_Multisim_expression_4():
    i = """task {
        task = "desmond:auto"
    }
    simulate {
        title = "NPT and no restraints, 24ps"
        effect_if = [[ "@*.*.annealing" ] 'annealing = off temperature = "@*.*.temperature[0][0]"' ]
        time = 24
        ensemble = {
            class = NPT
            method = Langevin
            thermostat.tau = 0.1
            barostat.tau = 2.0
        }
        eneseq.interval = 0.3
        trajectory.center = solute
    }
    simulate {
        meta = {
            cv = [
                { 
                    atom = [0 0 0 0] 
                    type = dihedral 
                    wall = 0  
                    width = 5.0 
                }
            ]
            cv_name = "$JOBNAME$[_replica$REPLICA$].cvseq"
            first = 0.0
            height = 0.03
            interval = 0.09
            name = "$JOBNAME$[_replica$REPLICA$].kerseq"
        }
    }
    """

    D = Multisim(string=i)
    
    assert len(D.dot.simulate[0].effect_if) == 2
    assert len(D.dot.simulate[0].effect_if[0]) == 1
    assert D.dot.simulate[0].ensemble.method == 'Langevin'
    assert D.dot.simulate[0].ensemble.barostat.tau == str(2.0)
    assert D.dot.simulate[0].eneseq.interval == str(0.3)
    assert D.dot.simulate[0].trajectory.center == 'solute' # <---
    
    assert len(D.dot.simulate[1].meta.cv) == 1
    assert len(D.dot.simulate[1].meta.cv[0].atom) == 4
    assert D.dot.simulate[1].meta.height == str(0.03)
    assert D.dot.simulate[-1].meta.cv[0].wall == str(0)
    
    D.dot.simulate[0].ensemble.barostat.tau = 2.1
    assert D.dot.simulate[0].ensemble.barostat.tau == 2.1

    D.dot.simulate[1].meta.height = 0.06
    assert D.dot.simulate[1].meta.height == 0.06
    
    D.dot.simulate[-1].meta.cv[0].wall = 50.0
    assert D.dot.simulate[-1].meta.cv[0].wall == 50.0
    

def test_Multisim_expression_5():
    D= Multisim(template="desmond-md.msj")
    D.write()