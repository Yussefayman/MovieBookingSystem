package moviebookingsystem;
import database.mysql;
import java.sql.ResultSet;
import java.sql.SQLException;


public class movie {
   String name;
    String release_date;
    String genre;
    String description;
    String poster;
    String city,cinema;
    int hall_id;
    public movie(){
        
    }
    public movie(String name, String release_date, String genre, String description,String poster)
    {
        this.name = name;
        this.release_date = release_date;
        this.genre = genre;
        this.description = description;
        this.poster=poster;
        
    }
    public void save_movie(String name_, String release_date_, String genre_, String description_,String poster) throws SQLException{
        boolean flag=false;
        ResultSet rs= mysql.RunOnQuery("select * from movie");
        while(rs.next()){
            
            name= rs.getString("movie_name");
             release_date= rs.getString("movie_release_date");
             genre= rs.getString("movie_genre");
             description = rs.getString("movie_desc");
             int movie_id= Integer.valueOf(rs.getString("movie_id"));
             if(name.equals(name_)&&release_date.equals(release_date_)&&genre.equals(genre_)&&description.equals(description_)){
                 flag=true;
             }
        }
        if(flag==true){
        }
        else{
           mysql.runnonquery("insert into movie (movie_name,movie_release_date,movie_genre,movie_desc,movie_poster) values ('"+name_+"','"+release_date_+"','"+genre_+"','"+description_+"','"+poster+"')");
                  
      }
      
    }
    public int movie_id(String name) throws SQLException{
        ResultSet rs= mysql.RunOnQuery("select * from movie where movie_name='"+name+"'");
        while(rs.next()){
            String movie_id= rs.getString("movie_id");
            return Integer.valueOf(movie_id);
           
        }
      return 0;
    }
}

