from django.shortcuts import render, redirect
from django.contrib import messages  # import 'messages' function to display different messages to user upon form completion
from .forms import UserRegisterForm

def register(request):
    # if request submits a form with POST data, create a form with the data from the POST
    if request.method == 'POST':
        form = UserRegisterForm(request.POST)
        if form.is_valid():
            form.save()
            username = form.cleaned_data.get('username')
            messages.success(request, f'Account created for {username}!')
            return redirect('blog-home')
    # if request has no POST data, create blank form
    else:
        form = UserRegisterForm()

    return render(request,'users/register.html',{'form': form})