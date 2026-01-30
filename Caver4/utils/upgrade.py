



import json
import urllib
from urllib.error import HTTPError, URLError

from ..caver_pymol import ROOT_LOGGER,VERSION

logging=ROOT_LOGGER.getChild('upgrade')

def get_github_repo_tags(repo_url) -> list[str]:
    """
    Retrieve all released tags of a GitHub repository using urllib.

    Usage:
        tags = get_github_repo_tags("https://github.com/BradyAJohnston/MolecularNodes")
        print(tags)

    Args:
        repo_url (str): The URL of the GitHub repository.

    Returns:
        list: A list of tag names for the repository.
    """

    # Extract the owner and repo name from the URL
    parts = repo_url.split("/")
    owner = parts[-2]
    repo = parts[-1]

    # GitHub API URL for listing tags
    api_url = f"https://api.github.com/repos/{owner}/{repo}/tags"

    try:
        # Send a GET request to the GitHub API
        with urllib.request.urlopen(api_url) as response:
            # Read the response and decode from bytes to string
            response_data = response.read().decode()
            # Parse JSON response data
            tags = json.loads(response_data)
            # Extract the name of each tag
            tag_names = [tag["name"] for tag in tags]
            return tag_names
    except HTTPError as e:
        # Handle HTTP errors (e.g., repository not found, rate limit exceeded)
        logging.warning(f"GitHub API returned status code {e.code}")
        return []
    except URLError as e:
        # Handle URL errors (e.g., network issues)
        logging.error(f"Failed to reach the server. Reason: {e.reason}")
        return []
    

def has_updates() -> bool:
    url='https://github.com/YaoYinYing/caver-pymol-plugin'
    tags=get_github_repo_tags(url)

    if tags:
        latest_tag = tags[0]
        return latest_tag.lstrip('v') > VERSION
    else:
        return False
    
