from django import forms
import re


def is_cg_code(cg_code):
    return bool(re.match(r'cg[0-9]{8}', cg_code))


class InputForm(forms.Form):
    cg_code = forms.CharField(label='', min_length=10, max_length=10)

    def clean(self):
        cleaned_data = super().clean()
        cg_code = cleaned_data.get('cg_code')
        cg_code = is_cg_code(cg_code)

        if not cg_code:
            raise forms.ValidationError('Invalid IlmnID')
