from django.shortcuts import render, redirect
from django.contrib import messages  # import 'messages' function to display different messages to user upon form completion
from .forms import UserRegisterForm
from django.contrib.auth import logout
from django.contrib.auth.decorators import login_required


def logout_view(request):
    logout(request)
    return render(request, 'users/logout.html')


def register(request):
    # if request submits a form with POST data, create a form with the data from the POST
    if request.method == 'POST':
        form = UserRegisterForm(request.POST)
        if form.is_valid():
            form.save()
            username = form.cleaned_data.get('username')
            messages.success(request, f'Your account has been created! You are now able to log in!')
            return redirect('login')
    # if request has no POST data, create blank form
    else:
        form = UserRegisterForm()

    return render(request,'users/register.html',{'form': form})

@login_required  # a decorator adds functionality to functions
def profile(request):
    return render(request, 'users/profile.html')