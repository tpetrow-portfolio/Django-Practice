from django.shortcuts import render

posts = [
    {
        'author': 'Tyler',
        'title': 'Blog Post 1',
        'content': 'First Post Content',
        'date_posted': 'September 18, 2024'
    },
    {
        'author': 'Taylor',
        'title': 'Blog Post 2',
        'content': 'Second Post Content',
        'date_posted': 'September 19, 2024'
    },
]

def home(request):
    context = {  # context is a dictionary, the keys are accessable in our template
        'posts': posts
    }
    return render(request, 'blog/home.html', context)  # adding 'context' to render allows the rendered webpage to access the content dictionary

def about(request):
    return render(request, 'blog/about.html', {'title': 'About'})  # add a webpage title called "About" 