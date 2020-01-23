from findley_evaluation_functions import evaluate_findley


class Critera:
    def __init__(self, name, evaluation_function):
        self.name = name
        self.evaluation_function = evaluation_function


Findley = Critera(name='Findley', evaluation_function=evaluate_findley)
