from django.shortcuts import render, redirect
from django.contrib import messages  # import 'messages' function to display different messages to user upon form completion
from .forms import UserRegisterForm, UserUpdateForm, ProfileUpdateForm
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
    if request.method == 'POST':
        u_form = UserUpdateForm(request.POST, instance=request.user)
        p_form = ProfileUpdateForm(request.POST, request.FILES, instance=request.user.profile)
        
        if u_form.is_valid() and p_form.is_valid():
            u_form.save()
            p_form.save()
            messages.success(request, f'Your account has been updated!')
            return redirect('profile')
        
    else:
        u_form = UserUpdateForm(instance=request.user)
        p_form = ProfileUpdateForm(instance=request.user.profile)

    context = {
        'u_form': u_form,
        'p_form': p_form
    }

    return render(request, 'users/profile.html', context)