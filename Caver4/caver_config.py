# CAVER Copyright Notice
# ============================
#

'''
"THE BEERWARE LICENSE" (Revision 42):

Yinying wrote the refactored code.
As long as you retain this notice, you can do whatever you want with this stuff.
If we meet someday, and you think this stuff is worth it, you can buy me a beer
in return.
-- Yinying Yao

'''

from dataclasses import dataclass
import json
import logging
import os
from typing import Any, Optional

THIS_DIR = os.path.dirname(__file__)
CONFIG_TXT = os.path.join(THIS_DIR, "config", "config.txt")

@dataclass
class CaverConfig:

    # connect to main ui widgets
    output_dir: str = ""
    selection_name: str = ''

    start_point_x: float = 0.0
    start_point_y: float = 0.0
    start_point_z: float = 0.0


    # connect to config widgets
    customized_java_heap: int = 6000

    shell_radius: float = 3.0
    shell_depth: int = 2
    probe_radius: float = 0.7
    clustering_threshold: float = 1.5

    number_of_approximating_balls: int = 4
    
    max_distance: float = 4.0
    desired_radius: float = 1.8

    def has(self, key: str) -> bool:
        return hasattr(self, key)

    def delete(self, key: str):
        delattr(self, key)

    def get(self, key: str) -> Any:
        return getattr(self, key)

    def set(self, key: str, value: Any):
        setattr(self, key, value)

    @classmethod
    def from_json(cls, json_file: str) -> 'CaverConfig':
        """
        Load configuration information from a JSON file and create a CaverConfig object.
    
        This method is a class method, which means it can be called directly on the class rather than an instance of the class.
        It takes a path to a JSON file as input and returns an instance of CaverConfig initialized with the data from the JSON file.
    
        Parameters:
        - json_file: str - The path to the JSON file containing the configuration information.
    
        Returns:
        - CaverConfig: An instance of the CaverConfig class initialized with the data from the JSON file.
        """
        # Load JSON data from the specified file path
        json_data = json.load(open(json_file))
        # filter dataclass keys with fixed typing
        # `__match_args__` has all dataclass keys
        # `__dataclass_fields__` has all dataclass key fields
        # Initialize a new instance of the class using the data from the JSON, converting types as necessary
        new_self = cls(**{k: cls.__dict__['__dataclass_fields__'][k].type(v)
                       for k, v in json_data.items() if k in cls.__dict__['__match_args__']})
        # add extra keys that come from json
        # Iterate through the JSON data again, adding any extra keys not included in the initial dataclass keys
        for k, v in json_data.items():
            if k in cls.__dict__['__match_args__']:
                continue
            new_self.set(k, v)
        return new_self

    def to_json(self, json_file: str):
        json.dump(self.__dict__, open(json_file, 'w'))

    def set_value(self, key: str, value: Any):
        setattr(self, key, value)

    @classmethod
    def from_txt(cls, txt_file: Optional[str] = None):
        '''
        Load configuration from a specified text file or the default configuration file if none is provided.
        
        This method reads the configuration from a text file, where each line represents a configuration item.
        The format of each line is "key value", with the following rules:
        - Lines starting with "#" are comments and are ignored.
        - If "#" appears in the middle of a line, the part after "#" is treated as a comment and is ignored.
        - Boolean values: "yes" represents True, "no" represents False.
        - String values: If they contain spaces, they must be enclosed in double quotes "".
        - List values: Space-separated.
        
        Parameters:
        - txt_file (Optional[str]): The path to the text file containing the configuration. If not provided, the default configuration file is used.
        
        Returns:
        - An instance of the class initialized with the configuration read from the file.
        '''
        # Initialize a new instance of the class
        new_self = cls()
        
        # Open the specified configuration file or the default configuration file
        with open(txt_file or CONFIG_TXT) as f:
            # Read each line of the file
            for l in f.readlines():
                # Remove leading and trailing whitespace
                l = l.strip()
                
                # If the line contains "#" and does not start with "#", remove the comment part
                if '#' in l and not l.startswith('#'):
                    l = l[0:l.rfind("#") - 1]
                    # remove trailing whitespaces
                    l = l.rstrip(' ')
                
                # If the line starts with "#" or is empty, skip it
                if l.startswith("#") or not l:
                    continue
                
                # Split the line into key and value parts
                parsed = l.split(' ', 1)
                key = parsed[0]
                # If the key part is empty or there is no value part, skip it
                if len(parsed) <= 1:
                    logging.debug('skipping ' + key)
                    continue
                
                # Set the value for the key in the new instance
                new_self._set_value(key, parsed[1])
        
        # Return the initialized new instance
        return new_self

    def _set_value(self, key: str, new_value: str):
        if hasattr(self, key):
            v_type = type(getattr(self, key))
            # align to the self variable type
            setattr(self, key, v_type(new_value) if v_type != bool else CaverConfig._true_or_false(new_value))
            return
        # unknown variable, just set it as str
        setattr(self, key, new_value)

    def _get_value(self, key: str) -> str:
        # if the variable is a defined slot in this dataclass, return with the correct type
        value = getattr(self, key)
        v_type = type(value)

        # if it's a boolean, return yes or no according to the config.txt
        if v_type == bool:
            return CaverConfig._yes_or_no(value)
        return str(value)

    @staticmethod
    def _yes_or_no(v: bool) -> str:
        return 'yes' if v else 'no'

    @staticmethod
    def _true_or_false(v: str) -> bool:
        return True if v == 'yes' else False

    def to_txt(self, txt_file: str):
        '''
        Export the current configuration to a TXT file for the JAR file to read.
        The format of the TXT file is as follows:
        - Boolean values are represented by "True" or "False"
        - List values are separated by spaces
    
        Parameters:
        - txt_file (str): The path to the TXT file to be exported
    
        This function does not return any value.
        '''
        # Set the starting point coordinates if any of the start_point_x, start_point_y, or start_point_z attributes are not 0
        if any(getattr(self, f'start_point_{a}') != 0 for a in 'xyz'):
            self.set('starting_point_coordinates', " ".join(str(getattr(self, f'start_point_{a}')) for a in 'xyz'))
            logging.debug(f"starting point: {self.get('starting_point_coordinates')}")
    
        # Initialize a new list to store the updated configuration content
        new_txt_contents = []
        # Read the configuration template file
        with open(CONFIG_TXT) as template_f:
            template = template_f.readlines()
        for l in template:
            l = l.strip()
            
            # skip the template line note
            if 'used as a template' in l:
                continue
            # Directly append comment lines and empty lines to the new configuration content
            if l.startswith('#') or not l:
                new_txt_contents.append(l)
                continue
            # Remove everything after the last occurrence of the # character and trailing whitespaces for lines containing inline comments
            if '#' in l and not l.startswith('#'):
                l = l[0:l.rfind("#") - 1]
                l = l.rstrip(' ')
    
            # Split each line into a key-value pair
            parsed = l.split(' ', 1)
            key = parsed[0]
            # Skip lines with no value set
            if len(parsed) <= 1:
                logging.debug('skipping ' + key)
            # Update keys that exist in the TXT file but not in the current configuration
            elif hasattr(self, key):
                self_val = getattr(self, key)
                _val_type = type(self_val)
                # Handle boolean values
                if _val_type == bool:
                    new_txt_contents.append(f'{key} {CaverConfig._yes_or_no(self_val)}')
                # Skip unset values
                elif self_val == '???':
                    continue
                # Handle float values
                elif _val_type == float:
                    new_txt_contents.append(f'{key} {self_val:.1f}')
                # Handle other types of values
                else:
                    new_txt_contents.append(f'{key} {str(self_val)}')
    
        # Create the directory for the TXT file if it does not exist
        os.makedirs(os.path.dirname(txt_file), exist_ok=True)
        # Write the updated configuration content to the TXT file
        with open(txt_file, 'w') as f:
            f.write('\n'.join(new_txt_contents))
