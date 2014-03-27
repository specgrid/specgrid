

class CompositeModel(object):

    def __init__(self, models_list):
        self.models_list = models_list
        self.param2model = {}
        for model in models_list:
            self.param2model.update(dict([(param, model) for param in model.parameters]))

    def __getattr__(self, item):
        if item in self.param2model:
            return getattr(self.param2model[item], item)
        else:
            raise AttributeError

    def __setattr__(self, item, value):
        if item in self.param2model:
            return setattr(self.param2model[item], item, value)
        else:
            raise AttributeError
