from django.shortcuts import render
from .models import Post


def home(request):
    context = {  # context is a dictionary, the keys are accessable in our template
        'posts': Post.objects.all()  # query the database
    }
    return render(request, 'blog/home.html', context)  # adding 'context' to render allows the rendered webpage to access the content dictionary

def about(request):
    return render(request, 'blog/about.html', {'title': 'About'})  # add a webpage title called "About" 