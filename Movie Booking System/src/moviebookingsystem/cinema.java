package moviebookingsystem;

import database.mysql;
import java.sql.ResultSet;
import java.sql.SQLException;


public class cinema {
    String name;
   String city;
   movie m = new movie();
   hall h= new hall();
    public cinema(String cinema, String city) {
        this.name=name;
        this.city=city;
    }
    public cinema(){
        
    }

   
    public void save_cinema(String cinema,String city) throws SQLException{
        boolean flag=false;
         ResultSet rs= mysql.RunOnQuery("select * from cinema");
       while(rs.next()){
           String cinema_id= rs.getString("cinema_id");
           String cinema_name= rs.getString("cinema_name");
           String cinema_city= rs.getString("cinema_city");
             
           if(cinema_name.equals(cinema)&&cinema_city.equals(city)){
             flag=true;
               
           }
          
       }
       if(flag==true){
           
       }
       else{
          
                mysql.runnonquery("insert into cinema (cinema_name,cinema_city) values ('"+cinema+"','"+city+"')");
       }
    }
    public boolean save_data(String cinema,String city,String movie_name, String movie_releasedate, String movie_genre, String movie_description,String movie_poster,int hallid) throws SQLException{
      
    

          int x=  cinema_id(cinema,city);
          if( h.check_hall(hallid,x))
          {
          
           m.save_movie(movie_name,  movie_releasedate, movie_genre, movie_description, movie_poster);
          int m_name= m.movie_id(movie_name);
         h.save_hall(hallid, x, movie_name);
          h.save_show(cinema, city, movie_name);
          
           return true;
          }
       
       return false;
    }
    public int cinema_id(String cinema_name,String cinema_city) throws SQLException{
        ResultSet rs= mysql.RunOnQuery("select * from cinema where cinema_name= '"+cinema_name+"' and cinema_city='"+cinema_city+"'");
       while(rs.next()){
           String cinema_id= rs.getString("cinema_id");
           return Integer.valueOf(cinema_id);
       }
       return 0;
    }
}
