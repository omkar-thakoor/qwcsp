#!/usr/bin/python3

# Generate random WCSP instances.

import random
import itertools

integer_cost = True
cost_range = (0, 20)

def generate_instance(num_variables, num_constraints, out_file):
    
    print('unknown {} 2 {} 99999'.format(num_variables, num_constraints), file=out_file)
    #print(*[random.randint(1,2) for _ in range(num_variables)], file=out_file)
    print('1 ' * num_variables, file=out_file)
    print('2 ' * num_variables, file=out_file)
    print('starting constraints of ' + str(num_constraints) + ' on ' + str(num_variables))

    constraint_set = set()

    for _ in range(num_constraints):
        while True:
            arity = random.randint(1, min(3, num_variables))
            variables = frozenset(random.sample(range(num_variables), arity))

            # Don't duplicate.
            if variables not in constraint_set:
                break

        print(arity, end=' ', file=out_file)

        constraint_set.add(variables)
        print(' '.join(str(i) for i in variables), end=' ', file=out_file)
        print('0 {}'.format(2 ** arity), file=out_file)

        for assignments in itertools.product(*([[0, 1],] * arity)):
            print(' '.join(str(i) for i in assignments), end=' ', file=out_file)
            if integer_cost:
                print(random.randint(*cost_range), file=out_file)
            else:
                print(random.uniform(*cost_range), file=out_file)

if __name__ == "__main__":
    num_instances = 4
    num_variables_range = (12,17)
    num_constraints_range = (20, 30)
    for inst in range(18, 18+num_instances):
        num_variables = random.randint(*num_variables_range)
        while True:
            num_constraints = random.randint(*num_constraints_range)
            if num_constraints < ((2 ** num_variables) * 0.5) :
                break
        out_file = open('instance' + str(inst) + '.txt', 'w+')
        generate_instance(num_variables, num_constraints, out_file)
        out_file.close()
