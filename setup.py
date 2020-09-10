from setuptools import setup


scripts=['NCResNet.py']

setup(
    name='NCResNet-wrapper',
    version='0.1.0',
    packages=['ncresnet'],
    scripts=scripts,
    include_package_data=True,
    package_data={'': ['data/*.*']},
    url='',
    license='',
    author='CHEN Yuelong',
    author_email='yuelong.chen.btr@gmail.com',
    description='NCResNet-wrapper'
)
