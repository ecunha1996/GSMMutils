import os.path
import time
from urllib.request import Request, urlopen
from bs4 import BeautifulSoup
import concurrent.futures

def test():
    from urllib.request import Request, urlopen
    from bs4 import BeautifulSoup
    import time
    use_cached = True
    # Base URL
    base_url = "https://www.futbin.com/players/evolutions"

    # Query Parameters
    params = (
        #"&position=CAM%2CRM%2CLM%2CST%2CRW%2CLW"
        "&player_rating=80-99"
        "&activeEvos=123%2C156%2C164%2C120"   #  rigid backbone
        #"&league=308"
        "&showSelectedEvos=true"
        "&sort=TotalIGS"
        "&order=desc"
        "&gender=men"
        #"&Pace=90"
    )
    rating_map = {}
    # Loop through multiple pages
    for page in range(1, 35):  # Adjust the range as needed
        url = f"{base_url}?page={page}{params}"  # Use `?page=<number>` followed by `params`
        print(url)
        req = Request(
            url,
            headers={"User-Agent": "Mozilla/5.0"}
        )

        try:
            if not os.path.exists('tmp'):
                os.mkdir('tmp')
            if not os.path.exists(f'webpage_{page}.html') or not use_cached:
                response = urlopen(req)
                html = response.read().decode("utf-8")  # Decode the HTML

                with open(f'tmp/webpage_{page}.html', 'w', encoding="utf-8") as f:
                    f.write(html)
            else:
                with open(f'tmp/webpage_{page}.html', 'r', encoding="utf-8") as f:
                    html = f.read()

            # Parse the HTML content
            soup = BeautifulSoup(html, "html.parser")

            # Extract player names
            player_names = [name.text.strip() for name in soup.select('.table-player-name')]

            # Extract player ratings
            player_ratings = [rating.text.strip() for rating in soup.select('.player-rating-card-text.font-standard.bold')]

            in_game_stats = [stats.text.strip() for stats in soup.select('.table-in-game-stats')][1:]

            evolutions = [evo.text.strip() for evo in soup.select('.evolutions-row')]

            # Print the results for the current page
            # print(f"Page {page}:")
            for name, rating, igs, evo in zip(player_names, player_ratings, in_game_stats, evolutions):
                # print(f"Player: {name}, Rating: {rating}, Stats: {igs}")
                rating_normalized = round(int(igs) / (int(rating)**1), 2)
                rating_map[name] = (rating_normalized, evo)

            time.sleep(1)

        except Exception as e:
            print(f"Error fetching page {page}: {e}")
            break
    # get top 10
    top_10 = sorted(rating_map.items(), key=lambda x: x[1], reverse=True)[0:10]
    #get all evos with Brazil in the name (the name is the second value in the tuple)
    brazil = {player: rating_map[player] for player in rating_map if 'Attack' in rating_map[player][1]}
    # sort by rating
    brazil = sorted(brazil.items(), key=lambda x: x[1], reverse=True)[0:10]
    print(brazil)

    print(top_10)

# def process_page(page, base_url, params, use_cached, attempt=0):
#     url = f"{base_url}?page={page}{params}"
#     print(url)
#     req = Request(
#         url,
#         headers={"User-Agent": "Mozilla/5.0"}
#     )
#
#     try:
#         if not os.path.exists('tmp'):
#             os.mkdir('tmp')
#         if not os.path.exists(f'tmp/webpage_{page}.html') or not use_cached:
#             response = urlopen(req)
#             html = response.read().decode("utf-8")  # Decode the HTML
#
#             with open(f'tmp/webpage_{page}.html', 'w', encoding="utf-8") as f:
#                 f.write(html)
#         else:
#             with open(f'tmp/webpage_{page}.html', 'r', encoding="utf-8") as f:
#                 html = f.read()
#
#         # Parse the HTML content
#         soup = BeautifulSoup(html, "html.parser")
#
#         # Extract player names
#         player_names = [name.text.strip() for name in soup.select('.table-player-name')]
#
#         # Extract player ratings
#         player_ratings = [rating.text.strip() for rating in soup.select('.player-rating-card-text.font-standard.bold')]
#
#         in_game_stats = [stats.text.strip() for stats in soup.select('.table-in-game-stats')][1:]
#
#         evolutions = [evo.text.strip() for evo in soup.select('.evolutions-row')]
#
#         rating_map = {}
#         for name, rating, igs, evo in zip(player_names, player_ratings, in_game_stats, evolutions):
#             rating_normalized = round(int(igs) / (int(rating)**1.5), 2)
#             rating_map[name] = (rating_normalized, evo)
#
#         return rating_map
#
#     except Exception as e:
#         if attempt < 10:
#             time.sleep(2)
#             return process_page(page, base_url, params, use_cached, attempt + 1)
#         else:
#             print(f"Error fetching page {page}: {e}")
#             return {}
#
# def test():
#     base_url = "https://www.futbin.com/players/evolutions"
#     params = (
#         "&position=CAM%2CRM%2CLM%2CST%2CRW%2CLW"
#         "&player_rating=40-90"
#         "&activeEvos=123%2C155%2C156"
#         "&showSelectedEvos=true"
#         "&sort=TotalIGS"
#         "&order=desc"
#         "&gender=men"
#         # "&Pace=90"
#     )
#     use_cached = False
#     rating_map = {}
#
#     with concurrent.futures.ThreadPoolExecutor() as executor:
#         futures = [executor.submit(process_page, page, base_url, params, use_cached) for page in range(1, 35)]
#         for future in concurrent.futures.as_completed(futures):
#             result = future.result()
#             rating_map.update(result)
#
#     top_10 = sorted(rating_map.items(), key=lambda x: x[1], reverse=True)[0:10]
#     brazil = {player: rating_map[player] for player in rating_map if 'Brazil' in rating_map[player][1]}
#     brazil = sorted(brazil.items(), key=lambda x: x[1], reverse=True)[0:10]
#     print(brazil)
#     print(top_10)



if __name__ == '__main__':
    # get_data()
    test()