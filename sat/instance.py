import re
from sortedcontainers import SortedSet


class SATInstance(object):
    def parse(self, line):
        line = re.sub(r'\s', '', line).lower()
        for clause in line.split('&'):
            literals_raw = clause.replace("(", "").replace(")", "").split("|")  # [x1,x2...]
            literals = {literal.replace("~", ""): 0 if literal.startswith("~") else 1 for literal in
                        literals_raw}
            if any(not literal.startswith('x') for literal in literals.keys()):
                raise ValueError("All literals must start with x!")
            self.variables.update(literals.keys())  # {x1, ... , xn}
            self.clauses.append(literals)

    def __init__(self):
        self.variables = SortedSet()
        self.clauses = []

    @property
    def number_of_variables(self):
        return len(self.variables)

    @property
    def number_of_clauses(self):
        return len(self.clauses)

    @classmethod
    def from_file(cls, file_path):
        instance = cls()
        with open(file_path) as file:
            for line in file.readlines():
                line = line.strip()
                if len(line) > 0 and not line.startswith('#'):
                    instance.parse(line)
        return instance

    @staticmethod
    def literal_to_string(literal, positive):
        return literal if positive else f"~{literal}"

    @staticmethod
    def to_assignment_graph_file(variables, file_path):
        graph_edges = []
        for i, var in enumerate(variables):
            literal, negated_literal = var, f'~{var}'
            a_start, a_end = f"a{str(i)}", f"a{str(i+1)}"
            graph_edges.extend([(a_start, literal), (a_start, negated_literal), (literal, a_end), (negated_literal, a_end)])

        with open(file_path, 'w+') as file:
            file.writelines("\n".join(" ".join(edge) for edge in graph_edges))
