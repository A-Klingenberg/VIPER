class _Singleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        """
        Return an instance of this class if it already exists, else construct one given the arguments. Also construct a
        new one when the kwarg 'force_refresh_singleton' is given.

        !!! Therefore, you cannot use the kwarg 'force_refresh_singleton' as part of the constructor of a subclass of this !!!

        :param args: Constructor args to pass along
        :param kwargs: Constructor kwargs to pass along
        :return: An instance of this class
        """
        force = kwargs.pop("force_refresh_singleton", False)
        if cls not in cls._instances or force:
            cls._instances[cls] = super(_Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]

    def get_instance(cls):
        """
        Returns an instance of this class, else raises an error.

        :return: An instance of this class
        """
        if c := cls._instances.get(cls, None):
            return c
        else:
            # Modify if you want to run custom logic if someone tries to get an instance, but none exist yet
            print(f"Trying to use {cls} before it has been initialized!")
            raise RuntimeError(f"Trying to use {cls} before it has been initialized!")
