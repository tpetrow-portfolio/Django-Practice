from django.db import models
from django.contrib.auth.models import User

class Profile(models.Model):
    # one to one relationship means one user for one profile, and one profile for one user
    # if we delete the profile, delete the user, but if we delete the user, don't delete the profile
    user = models.OneToOneField(User, on_delete=models.CASCADE)
    image = models.ImageField(default='default.jpg', upload_to='profile_pics')

    def __str__(self):
        return f'{self.user.username} Profile'

