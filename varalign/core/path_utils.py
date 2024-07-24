# varalign/core/path_utils.py

import os

def get_project_root():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))

def get_test_data_path():
    return os.path.join(get_project_root(), 'tests', 'data')

def get_issues_path():
    return os.path.join(get_project_root(), 'tests', 'issues')

if __name__ == '__main__':
    print(get_project_root())
    print(get_test_data_path())
    print(get_issues_path())
